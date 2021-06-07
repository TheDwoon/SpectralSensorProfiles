#include <iostream>
#include <utility>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include <Eigen/Dense>

#include "cie1931.h"
#include "rgb2spec.h"
#include "spectrum.h"
#include "xy_conversion.h"
#include "sampler.h"
#include "camera_xy_lookup.h"
#include "ppm_image.h"

#define C_PARAMS 6

// taken from rgb2spec.cpp of Sigmoid Paper
// rewritten to support ceres template for automatic derivatives
template<typename T, int N>
T eval_precise(const T* c, const T &lambda) {
    T x = c[0];
    for (int i = 1; i < N; i++) {
        x = x * lambda + c[i];
    }
    T y = 1.0 / ceres::sqrt(x * x + 1.0);
    return 0.5 * x * y + 0.5;
}

template<int N>
struct SpecIntegrator {
    SpecIntegrator(spec::Spectrum spec, Eigen::Vector3d target) : target_(std::move(target)), spectrum_(std::move(spec)) {}

    template<typename T>
    T integrate(const T* c) const {
        T integral(0);
        for (int i = 0; i < spectrum_.size; ++i) {
            T lambda(static_cast<double>(i) / (spectrum_.size - 1));
            T eval = eval_precise<T, N>(c, lambda);
            integral += spectrum_.spec[i] * eval;
        }

        return integral;
    }

    template<typename T>
    bool operator()(const T* parameter, T* residual) const {
        const T* r = parameter + N * 0;
        const T* g = parameter + N * 1;
        const T* b = parameter + N * 2;

        T red = integrate(r);
        T green = integrate(g);
        T blue = integrate(b);

        residual[0] = T(target_(0)) - red / (red + green + blue);
        residual[1] = T(target_(1)) - green / (red + green + blue);
        return true;
    }

    private:
        const Eigen::Vector3d target_;
        spec::Spectrum spectrum_;
};

#define ITERATIONS 1000
#define BATCHES 1500
#define EXAMPLES 50

#define REJECTION_DISTANCE 0.025
#define VECTOR_SAMPLES 200

void printVector(double* v, int n) {
    std::cout << "(";
    for (int i = 0; i < n; ++i) {
        std::cout << v[i];
        if (i < n - 1)
            std::cout << ", ";
    }

    std::cout << ")" << std::endl;
}

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);

    std::cout << "Loading xy converter... " << std::endl;
    spec::xy_conversion xyConverter = spec::loadConversion();

    std::cout << "Reading matrix from std::cin" << std::endl;
    double matrix[9];
    for (int i = 0; i < 9; ++i)
        std::cin >> matrix[i];

    Eigen::Matrix3d xyz2cam;
    xyz2cam << matrix[0],matrix[1],matrix[2],matrix[3],matrix[4],matrix[5],matrix[6],matrix[7],matrix[8];
    xyz2cam /= 10000.0;

    std::cout << "xyz2cam\n";
    std::cout << xyz2cam << std::endl;

    Eigen::Matrix3d cam2xyz = xyz2cam.inverse();
    std::cout << "cam2xyz\n";
    std::cout << cam2xyz << std::endl;

    double cs[3 * C_PARAMS] { 0 };

    cs[0 * C_PARAMS + C_PARAMS - 3] = -3;
    cs[1 * C_PARAMS + C_PARAMS - 3] = -3;
    cs[2 * C_PARAMS + C_PARAMS - 3] = -3;

    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> rgb2xyz((double*)srgb_to_xyz);

    std::ofstream fitting_file("../cfa_fitting.dat", std::ofstream::out | std::ofstream::trunc);
    fitting_file << "# Batch\tIterations\tInitial Error\tFinal Error\n";

    std::ofstream sampling_file("../cfa_sampling.dat", std::ofstream::out | std::ofstream::trunc);
    sampling_file << "# x\ty\n";

    Sampler sampler;
    for (int batch = 0; batch < BATCHES; ++batch) {
        std::cout << "*************** BATCH " << batch << " ***************" << std::endl;
        size_t generatedSamples = 0;

        ceres::Problem problem;
        while (generatedSamples < EXAMPLES) {
            double d1 = sampler.sample(1. / 3., 0.05, 3);
            double d2 = sampler.sample(1. / 3., 0.05, 3);
            if (d1 + d2 > 1)
                continue;

            Eigen::Vector3d xyzColor(d1, d2, 1 - d1 - d2);
            xyzColor /= xyzColor.sum();
            Eigen::Vector3d camColor = xyz2cam * xyzColor;
            camColor /= camColor.sum();

            spec::coeff_t nc = xyConverter.access(xyzColor(0), xyzColor(1));
            if (nc[0] == 0. && nc[1] == 0. && nc[2] == 0.)
                continue;

            sampling_file << xyzColor(0) << "\t" << xyzColor(1) << "\n";

            generatedSamples += 1;
            spec::Spectrum spectrum(380, 780, nc.data(), 3);
            ceres::CostFunction *ct_fun = new ceres::AutoDiffCostFunction<SpecIntegrator<C_PARAMS>, 2, 3 * C_PARAMS>(
                    new SpecIntegrator<C_PARAMS>(spectrum, camColor));
            problem.AddResidualBlock(ct_fun, nullptr, cs);
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < C_PARAMS; ++j) {
                problem.SetParameterLowerBound(cs, i * C_PARAMS + j, -50.);
                problem.SetParameterUpperBound(cs, i * C_PARAMS + j, 50.);
            }
        }

        ceres::Solver::Options options;
        options.linear_solver_type = ceres::CGNR;
//        options.minimizer_progress_to_stdout = true;
        options.max_num_consecutive_invalid_steps = 5;
        options.max_num_iterations = ITERATIONS;
        options.num_threads = 12;
        options.gradient_tolerance = 1.e-12;
        options.function_tolerance = 1.e-12;
        options.parameter_tolerance = 1.e-12;

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << std::endl;

        // ignore in first batch. It's just ruining the scaling
        if (batch > 0) {
            fitting_file << batch << "\t";
            fitting_file << summary.iterations.size() << "\t";
            fitting_file << summary.initial_cost << "\t";
            fitting_file << summary.final_cost << "\n";
        }
    }

    fitting_file.flush();
    fitting_file.close();

    sampling_file.flush();
    sampling_file.close();

    spec::Spectrum camRed = spec::fromNormalized(380, 780, cs + 0 * C_PARAMS, C_PARAMS);
    std::cout << "Red: " << '\n';
    camRed.print(3);

    spec::Spectrum camGreen = spec::fromNormalized(380, 780, cs + 1 * C_PARAMS, C_PARAMS);
    std::cout << "Green: " << '\n';
    camGreen.print(3);

    spec::Spectrum camBlue = spec::fromNormalized(380, 780, cs + 2 * C_PARAMS, C_PARAMS);
    std::cout << "Blue: " << '\n';
    camBlue.print(3);

    camRed.store("../camera_red.bin");
    camGreen.store("../camera_green.bin");
    camBlue.store("../camera_blue.bin");

    for (int i = 0; i < 3; ++i) {
        printVector(cs + i * C_PARAMS, C_PARAMS);
    }

    // writing to output file
    assert(camRed.size == camGreen.size && camGreen.size == camBlue.size);
    assert(camRed.start == camGreen.start && camGreen.start == camBlue.start);

    std::ofstream spectrum_file("../cfa_spectra.dat", std::ofstream::out | std::ofstream::trunc);
    spectrum_file << "# Parameters for CFA in order red, green, blue\n";
    for (int i = 0; i < 3; ++i) {
        spectrum_file << "# " << cs[i * C_PARAMS + 0];
        for (int j = 1; j < C_PARAMS; ++j) {
            spectrum_file << "\t" << cs[i * C_PARAMS + j];
        }
        spectrum_file << "\n";
    }

    spectrum_file << "# Lambda\tRed\tGreen\tBlue\n";
    for (size_t i = 0; i < camRed.size; ++i) {
        spectrum_file << (i + camRed.start) << "\t" << camRed.spec[i] << "\t" << camGreen.spec[i] << "\t" << camBlue.spec[i] << "\n";
    }

    spectrum_file.flush();
    spectrum_file.close();

    std::ofstream vector_file("../cfa_vector.dat", std::ofstream::out | std::ofstream::trunc);
    vector_file << "# sourceX sourceY targetX targetY\n";

    // writing vector field data
    std::vector<Eigen::Vector3d> sampledPoints;

    int attempts = 10000;
    int vectorSamples = VECTOR_SAMPLES;
    while (vectorSamples > 0 && --attempts > 0) {
        const double d1 = sampler.sample(1. / 3., 0.33, 3);
        const double d2 = sampler.sample(1. / 3., 0.33, 3);
        if (d1 + d2 > 1.)
            continue;

        Eigen::Vector3d xyzColor = Eigen::Vector3d(d1, d2, 1 - d1 - d2);
        xyzColor /= xyzColor.sum();

        auto testFn = [&xyzColor](Eigen::Vector3d& v) {
            double dx = xyzColor(0) - v(0);
            double dy = xyzColor(1) - v(1);
            return dx * dx + dy * dy < REJECTION_DISTANCE * REJECTION_DISTANCE;
        };
        if (std::any_of(sampledPoints.begin(), sampledPoints.end(), testFn))
            continue;
        else
            sampledPoints.push_back(xyzColor);

        spec::coeff_t nc = xyConverter.access(xyzColor(0), xyzColor(1));
        if (nc[0] == 0. && nc[1] == 0. && nc[2] == 0.)
            continue;

        // cam color when integrating using fitted CFA
        spec::Spectrum spectrum(380, 780, nc.data(), 3);
        Eigen::Vector3d integratedColor(spec::integrateSpectrum(spectrum, camRed),
                                        spec::integrateSpectrum(spectrum, camGreen),
                                        spec::integrateSpectrum(spectrum, camBlue));

//        integratedColor /= integratedColor.sum();
        Eigen::Vector3d xyzIntegrated = cam2xyz * integratedColor;
        xyzIntegrated /= xyzIntegrated.sum();

        vector_file << xyzColor(0) << " " << xyzColor(1) << " ";
        vector_file << xyzIntegrated(0) << " " << xyzIntegrated(1) << "\n";

        vectorSamples -= 1;
    }

    vector_file.flush();
    vector_file.close();

    // Generate inverted function lookup
    const uint16_t lookupPrecision = 512;
    spec::camera_xy_lookup lookup = spec::camera_xy_lookup::generateXYLookUp(camRed, camGreen, camBlue, xyConverter, lookupPrecision);
    lookup.store();
    lookup.storeVectorFile(20, xyConverter, camRed, camGreen, camBlue);

//    if (argc >= 2) {
        std::string path("/home/daniel/git/masterarbeit/sensor_spectra/cameras/Kodak DCS460D_measured.dat");
//        std::string path("/home/daniel/git/masterarbeit/sensor_spectra/cameras/Canon EOS 5D Mark II_measured.dat");
        spec::Spectrum redRef(0, 0, { 0. });
        spec::Spectrum greenRef(0, 0, { 0. });
        spec::Spectrum blueRef(0, 0, { 0. });
        spec::parseFromFile(&redRef, &greenRef, &blueRef, path);
        lookup.storeVectorFile(20, xyConverter, redRef, greenRef, blueRef, "../calc_error_ref.dat");

        spec::camera_xy_lookup alternateLookup = spec::camera_xy_lookup::generateXYLookUp(redRef, greenRef, blueRef, xyConverter, lookupPrecision);
        alternateLookup.store("../camera_lookup_ref.bin");
//    }

    const double step = 1. / (lookupPrecision - 1);

    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> xyz2rgb((double*)xyz_to_srgb);

    // camera transformed rg_srgb_image
    const uint32_t rg_image_res = lookupPrecision;
    const uint32_t xy_image_res = lookupPrecision;

    // Contains an image with [rg] as coordinates and xyz as pixel color
    ppm_image rg_xyz_image(rg_image_res, rg_image_res);
    // contains an image with [rg] as coordinates and sRGB as pixel color
    ppm_image rg_srgb_image(rg_image_res, rg_image_res);
    // contains an image with [xy] as coordinates and camrgb as pixel color
    ppm_image xy_camrgb_image(xy_image_res, xy_image_res);

    for (int yIdx = 0; yIdx < lookupPrecision; ++yIdx) {
        const double y = yIdx * step;
        for (int xIdx = 0; xIdx < lookupPrecision; ++xIdx) {
            const double x = xIdx * step;
            auto [r, g] = lookup.accessRaw(x, y);
            if (r < 0 || g < 0)
                continue;

            Eigen::Vector3d xyz { x, y, 1. - x - y };
            Eigen::Vector3d rgb = xyz2rgb * xyz;
            rgb(0) = std::min(std::max(0., rgb(0)), 1.);
            rgb(1) = std::min(std::max(0., rgb(1)), 1.);
            rgb(2) = std::min(std::max(0., rgb(2)), 1.);
            if (rgb.sum() > 0)
                rgb /= rgb.sum();

            rg_srgb_image.write(static_cast<uint32_t>(r * rg_image_res), rg_image_res - static_cast<uint32_t>(g * rg_image_res),
                                static_cast<uint8_t >(rgb(0) * 255), static_cast<uint8_t>(rgb(1) * 255), static_cast<uint8_t>(rgb(2) * 255));
            rg_xyz_image.write(static_cast<uint32_t>(r * rg_image_res), rg_image_res - static_cast<uint32_t>(g * rg_image_res),
                               static_cast<uint8_t >(xyz(0) * 255), static_cast<uint8_t>(xyz(1) * 255), static_cast<uint8_t>(xyz(2) * 255));
            xy_camrgb_image.write(static_cast<uint32_t>(x * xy_image_res), xy_image_res - static_cast<uint32_t>(y * xy_image_res),
                                  static_cast<uint8_t >(r * 255.), static_cast<uint8_t>(g * 255.), static_cast<uint8_t>(255 - (r + g) * 255));
        }
    }

    rg_srgb_image.store("../rg_srgb.ppm");
    rg_xyz_image.store("../rg_xyz.ppm");
    xy_camrgb_image.store("../xy_camrgb.ppm");

    // contains an image with [rg] as coordinates and xyz as pixel color (function inverse!)
    ppm_image rg_inverse(lookupPrecision, lookupPrecision);
    for (int yIdx = 0; yIdx < lookupPrecision; ++yIdx) {
        const double g = yIdx * step;
        for (int xIdx = 0; xIdx < lookupPrecision; ++xIdx) {
            const double r = xIdx * step;

            auto [x, y] = lookup.access(r, g);
            if (x < 0. || y < 0.)
                rg_inverse.write(xIdx, lookupPrecision - 1 - yIdx, 0, 0, 0);
            else
                rg_inverse.write(xIdx, lookupPrecision - 1 - yIdx,
                                static_cast<uint8_t>(x * 255.), static_cast<uint8_t>(y * 255.), static_cast<uint8_t>(255. - (x + y) * 255.));
        }
    }

    rg_inverse.store("../rg_invserse.ppm");

    return 0;
}
