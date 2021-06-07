#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <numeric>
#include "camera_xy_lookup.h"
#include "spectrum.h"
#include "ppm_image.h"
#include "spectral_image.h"
#include "cie1931.h"
#include "color_correction.h"

/*
 * FIXME: This whole code here is cruel to use and just super hacky.
 * Rewrite this to use more methods and be easier to use.
 * Integrate spectral image into 3 double channels and use them to maximize stuff.
 * Put these util functions into a separate cpp file.
 */

const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> xyz2srgb((double*)xyz_to_srgb);
const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> xyz2rec2020((double*)xyz_to_rec2020);

struct image_d {
    uint32_t width;
    uint32_t height;

    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;

    void resize(uint32_t w, uint32_t h) {
        width = w;
        height = h;

        uint32_t size = w * h;
        a.resize(size);
        b.resize(size);
        c.resize(size);
    }
};

void saveErrorImage(uint32_t width, uint32_t height, double max_error, const std::vector<double> &delta_es, const std::string &path) {
    ppm_image image(width, height);

    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            double e = delta_es[y * width + x] / max_error;
            image.write(x, y, 255. * e, 0, 0);
        }
    }

    image.store(path);
}

Eigen::Vector3d applyColorSpaceTransform(const Eigen::Vector3d &color, const Eigen::Matrix3d &mat, double brightness) {
    Eigen::Vector3d result = mat * color;
    result /= result.sum();

    result *= brightness;
    result(0) = std::clamp(result(0), 0., 1.);
    result(1) = std::clamp(result(1), 0., 1.);
    result(2) = std::clamp(result(2), 0., 1.);
    result(0) = GAMMA(result(0));
    result(1) = GAMMA(result(1));
    result(2) = GAMMA(result(2));

    return result;
}

void writeImage(ppm_image &image, uint32_t x, uint32_t y, const Eigen::Vector3d color) {
    image.write(x, y, 255. * color(0), 255. *  color(1), 255. * color(2));
}

void printDeltaEInfo(const std::vector<double>& delta_es) {
    double min_delta_e = *std::min_element(delta_es.cbegin(), delta_es.cend());
    double max_delta_e = *std::max_element(delta_es.cbegin(), delta_es.cend());
    double avg_delta_e = std::accumulate(delta_es.cbegin(), delta_es.cend(), 0.f) / delta_es.size();

    std::cout << "min \u0394e " << min_delta_e << " max \u0394e " << max_delta_e << " avg \u0394e " << avg_delta_e << "\n";
}

void printColorDeltaInfo(const std::vector<double>& delta_es) {
    double min_delta_e = *std::min_element(delta_es.cbegin(), delta_es.cend());
    double max_delta_e = *std::max_element(delta_es.cbegin(), delta_es.cend());
    double avg_delta_e = std::accumulate(delta_es.cbegin(), delta_es.cend(), 0.f) / delta_es.size();

    std::cout << "min \u0394xy " << min_delta_e << " max \u0394xy " << max_delta_e << " avg \u0394xy " << avg_delta_e << "\n";
}

void sanitizeResponse(std::vector<double> &channel_a, std::vector<double> &channel_b, std::vector<double> &channel_c) {
    double max_a = *std::max_element(channel_a.begin(), channel_a.end());
    double max_b = *std::max_element(channel_b.begin(), channel_b.end());
    double max_c = *std::max_element(channel_c.begin(), channel_c.end());
    double max_response = std::max(max_a, std::max(max_b, max_c));

    assert(channel_a.size() == channel_b.size() && channel_b.size() == channel_c.size());
    for (size_t i = 0; i < channel_a.size(); ++i) {
        channel_a[i] /= max_response;
        channel_b[i] /= max_response;
        channel_c[i] /= max_response;
    }
}

image_d renderImage(const spec::Spectrum& red, const spec::Spectrum& green, const spec::Spectrum &blue, const spectral_image &spec_image) {
    const uint32_t size = spec_image.m_width * spec_image.m_height;

    // create resulting image with correct size
    image_d rendered_image;
    rendered_image.width = spec_image.m_width;
    rendered_image.height = spec_image.m_height;
    std::vector<double> &red_channel = rendered_image.a;
    std::vector<double> &green_channel = rendered_image.b;
    std::vector<double> &blue_channel = rendered_image.c;

    red_channel.resize(size);
    green_channel.resize(size);
    blue_channel.resize(size);

    double* buffer = new double[spec_image.m_dimension];
    for (uint32_t y = 0; y < spec_image.m_height; ++y) {
        for (uint32_t x = 0; x < spec_image.m_width; ++x) {
            const uint32_t idx = y * spec_image.m_width + x;
            spec_image.gather(buffer, x, y);
            spec::Spectrum pixel_spectrum(400, 700, 5, buffer);

            // integrate using provided spectrum
            red_channel[idx] = spec::interpolatedIntegration(red, pixel_spectrum);
            green_channel[idx] = spec::interpolatedIntegration(green, pixel_spectrum);
            blue_channel[idx] = spec::interpolatedIntegration(blue, pixel_spectrum);
        }
    }

    delete[] buffer;

    sanitizeResponse(red_channel, green_channel, blue_channel);
    return rendered_image;
}

std::vector<double> computeResponseBrightness(const image_d &image) {
    std::vector<double> brightness(image.width * image.height);

    for (uint32_t y = 0; y < image.height; ++y) {
        for (uint32_t x = 0; x < image.width; ++x) {
            const uint32_t idx = y * image.width + x;
            brightness[idx] = image.a[idx] + image.b[idx] + image.c[idx];
        }
    }

    return brightness;
}

void storeViewableImages(const image_d &xyz_image, const std::vector<double> &brightness, const std::string &path_xyz,
                       const std::string &path_srgb, const std::string &path_rec2020) {
    double max_brightness = *std::max_element(brightness.cbegin(), brightness.cend());

    ppm_image image_xyz(xyz_image.width, xyz_image.height);
    ppm_image image_srgb(xyz_image.width, xyz_image.height);
    ppm_image image_rec2020(xyz_image.width, xyz_image.height);

    for (uint32_t y = 0; y < xyz_image.height; ++y) {
        for (uint32_t x = 0; x < xyz_image.width; ++x) {
            const uint32_t idx = y * xyz_image.width + x;
            Eigen::Vector3d xyz { xyz_image.a[idx], xyz_image.b[idx], xyz_image.c[idx] };
            xyz /= xyz.sum();

            writeImage(image_xyz, x, y, xyz);

            Eigen::Vector3d srgb = applyColorSpaceTransform(xyz, xyz2srgb, brightness[idx] / max_brightness);
            writeImage(image_srgb, x, y, srgb);

            Eigen::Vector3d rec2020 = applyColorSpaceTransform(xyz, xyz2rec2020, brightness[idx] / max_brightness);
            writeImage(image_rec2020, x, y, rec2020);
        }
    }

    image_xyz.store(path_xyz);
    image_srgb.store(path_srgb);
    image_rec2020.store(path_rec2020);
}

image_d transformResponse_lookup(const image_d &response, const spec::camera_xy_lookup &lookup) {
    image_d xyz_image;
    xyz_image.resize(response.width, response.height);

    for (uint32_t y = 0; y < response.height; ++y) {
        for (uint32_t x = 0; x < response.width; ++x) {
            const uint32_t idx = y * response.width + x;

            double r = response.a[idx];
            double g = response.b[idx];
            double b = response.c[idx];
            double brightness = r + g + b;

            r /= brightness;
            g /= brightness;
            b /= brightness;

            std::array<double, 2> inverted = lookup.access(r, g);

            xyz_image.a[idx] = inverted[0];
            xyz_image.b[idx] = inverted[1];
            xyz_image.c[idx] = 1. - inverted[0] - inverted[1];
        }
    }

    return xyz_image;
}

image_d transformResponse_3x3(const image_d &response, const Eigen::Matrix3d &xyz2cam) {
    image_d xyz_image;
    xyz_image.resize(response.width, response.height);

    for (uint32_t y = 0; y < response.height; ++y) {
        for (uint32_t x = 0; x < response.width; ++x) {
            const uint32_t idx = y * response.width + x;

            Eigen::Vector3d cam_rgb { response.a[idx], response.b[idx], response.c[idx] };
            cam_rgb /= cam_rgb.sum();
            Eigen::Vector3d xyz = xyz2cam.inverse() * cam_rgb;
            xyz /= xyz.sum();

            xyz_image.a[idx] = xyz(0);
            xyz_image.b[idx] = xyz(1);
            xyz_image.c[idx] = xyz(2);
        }
    }

    return xyz_image;
}

std::vector<double> computeImageDeltaE(const image_d &a, const std::vector<double> &a_brightness,
                                  const image_d &b, const std::vector<double> &b_brightness) {
    assert(a.width == b.width);
    assert(a.height == b.height);

    double a_max_brightness = *std::max_element(a_brightness.cbegin(), a_brightness.cend());
    double b_max_brightness = *std::max_element(b_brightness.cbegin(), b_brightness.cend());

    std::vector<double> delta_es(a.width * a.height);

    for (uint32_t y_image = 0; y_image < a.height; ++y_image) {
        for (uint32_t x_image = 0; x_image < a.width; ++x_image) {
            const uint32_t idx = y_image * a.width + x_image;

            float a_xyz[3] { static_cast<float>(a.a[idx]), static_cast<float>(a.b[idx]), static_cast<float>(a.c[idx]) };
            float a_sum = a_xyz[0] + a_xyz[1] + a_xyz[2];
            a_xyz[0] = a_xyz[0] / a_sum * (a_brightness[idx] / a_max_brightness);
            a_xyz[1] = a_xyz[1] / a_sum * (a_brightness[idx] / a_max_brightness);
            a_xyz[2] = a_xyz[2] / a_sum * (a_brightness[idx] / a_max_brightness);

            float b_xyz[3] { static_cast<float>(b.a[idx]), static_cast<float>(b.b[idx]), static_cast<float>(b.c[idx]) };
            float b_sum = b_xyz[0] + b_xyz[1] + b_xyz[2];
            b_xyz[0] = b_xyz[0] / b_sum * (b_brightness[idx] / b_max_brightness);
            b_xyz[1] = b_xyz[1] / b_sum * (b_brightness[idx] / b_max_brightness);
            b_xyz[2] = b_xyz[2] / b_sum * (b_brightness[idx] / b_max_brightness);

            delta_es[idx] = compute_delta_e(a_xyz, b_xyz);
        }
    }

    return delta_es;
}

std::vector<double> computeImageColorDifference(const image_d &a, const image_d &b) {
    assert(a.width == b.width);
    assert(a.height == b.height);

    std::vector<double> deltas(a.width * a.height);
    for (uint32_t y_image = 0; y_image < a.height; ++y_image) {
        for (uint32_t x_image = 0; x_image < a.width; ++x_image) {
            const uint32_t idx = y_image * a.width + x_image;

            Eigen::Vector3d xyz_a { a.a[idx], a.b[idx], a.c[idx] };
            xyz_a /= xyz_a.sum();

            Eigen::Vector3d xyz_b { b.a[idx], b.b[idx], b.c[idx] };
            xyz_b /= xyz_b.sum();

            Eigen::Vector3d diff = xyz_b - xyz_a;

            // Color difference is measured in xy not xyz :)
            deltas[idx] = std::sqrt(diff(0) * diff(0) + diff(1) * diff(1));
        }
    }

    return deltas;
}

void visualizeChromaticError(const image_d &image, const std::vector<double> &delta_es, double max_error, const std::string &path) {
    const uint32_t resolution = 512;
    ppm_image res(resolution, resolution);

    for (uint32_t y = 0; y < image.height; ++y) {
        for (uint32_t x = 0; x < image.width; ++x) {
            const uint32_t idx = y * image.width + x;

            double relative_error = delta_es[idx] / max_error;
            double xyz[3] { image.a[idx], image.b[idx], image.c[idx] };
            double b = xyz[0] + xyz[1] + xyz[2];
            xyz[0] /= b;
            xyz[1] /= b;
            xyz[2] /= b;

            res.write(xyz[0] * resolution, resolution - 1. - (xyz[1] * resolution), 255. * relative_error, 0, 0);
        }
    }

    res.store(path);
}

int main(int argc, char** argv) {
    assert(argc >= 2);

    spec::Spectrum cam_red = spec::loadFromFile("../camera_red.bin");
    spec::Spectrum cam_green = spec::loadFromFile("../camera_green.bin");
    spec::Spectrum cam_blue = spec::loadFromFile("../camera_blue.bin");

    spec::Spectrum cie_x_spec = spec::loadFromFile("../cie1931_x.bin");
    spec::Spectrum cie_y_spec = spec::loadFromFile("../cie1931_y.bin");
    spec::Spectrum cie_z_spec = spec::loadFromFile("../cie1931_z.bin");

    spec::Spectrum red_ref(0, 0, { 0. });
    spec::Spectrum green_ref(0, 0, { 0. });
    spec::Spectrum blue_ref(0, 0, { 0. });

    if (argc == 3) {
        spec::parseFromFile(&red_ref, &green_ref, &blue_ref, argv[2]);
    }

    cam_red.print(3);
    cam_green.print(3);
    cam_blue.print(3);

    std::string path(argv[1]);
    std::cout << "Loading " << path << std::endl;
    spectral_image spec_image;
    spec_image.load(path);
    spec_image.gammaCorrection();

    std::cout << "Loaded spectral image: " << spec_image.m_width << "x" << spec_image.m_height << " [" << spec_image.m_dimension << "]" << std::endl;

    image_d cam_response = renderImage(cam_red, cam_green, cam_blue, spec_image);
    std::vector<double> cam_brightness_vector = computeResponseBrightness(cam_response);

    image_d cie_response = renderImage(cie_x_spec, cie_y_spec, cie_z_spec, spec_image);
    std::vector<double> cie_brightness_vector = computeResponseBrightness(cie_response);

    image_d ref_response = renderImage(red_ref, green_ref, blue_ref, spec_image);
    std::vector<double> ref_brightness_vector = computeResponseBrightness(ref_response);

    // render camrgb image
    ppm_image cam_image(spec_image.m_width, spec_image.m_height);
    for (uint32_t y = 0; y < spec_image.m_height; ++y) {
        for (uint32_t x = 0; x < spec_image.m_width; ++x) {
            const uint32_t idx = y * spec_image.m_width + x;
            cam_image.write(x, y, 255. * cam_response.a[idx], 255. * cam_response.b[idx], 255. * cam_response.c[idx]);
        }
    }

    cam_image.store("../spectral_camrgb.ppm");

    spec::camera_xy_lookup ref_lookup = spec::camera_xy_lookup::loadCameraXYLookUp("../camera_lookup_ref.bin");

    Eigen::Matrix3d xyz2cam;
    // Kodak DCS460D
    xyz2cam <<  1.0592, -0.2206, -0.0967, -0.1944, 1.1685, 0.023, 0.2206, 0.067,  0.1273;

    // Canon EOS 5D Mark II
//    xyz2cam << 4716., 603., -830., -7798., 15474., 2480., -1496., 1937., 6651.;
//    xyz2cam /= 10000.;

    spec::camera_xy_lookup lookup = spec::camera_xy_lookup::loadCameraXYLookUp();

    double max_error = 0.00001;

    // write CIE
    storeViewableImages(cie_response, cie_brightness_vector, "../cie_xyz.ppm", "../cie_srgb.ppm", "../cie_rec2020.ppm");

    // fitted spectra, fitted lookup
    image_d image_fitted_fitted = transformResponse_lookup(cam_response, lookup);
    storeViewableImages(image_fitted_fitted, cam_brightness_vector, "../spectra_xyz.ppm", "../spectra_srgb.ppm", "../spectra_rec2020.ppm");
    std::vector<double> cie_fitted_fitted_de = computeImageDeltaE(cie_response, cie_brightness_vector, image_fitted_fitted, cam_brightness_vector);
    max_error = std::max(max_error, *std::max_element(cie_fitted_fitted_de.cbegin(), cie_fitted_fitted_de.cend()));
    std::cout << "cie <> fitted_fitted: ";
    printDeltaEInfo(cie_fitted_fitted_de);

    // reference spectra, 3x3 matrix
    image_d image_ref_3x3 = transformResponse_3x3(ref_response, xyz2cam);
    storeViewableImages(image_ref_3x3, ref_brightness_vector, "../spectra_ref3x3_xyz.ppm", "../spectra_ref3x3_srgb.ppm", "../spectra_ref3x3_rec2020.ppm");
    std::vector<double> cie_ref3x3_de = computeImageDeltaE(cie_response, cie_brightness_vector, image_ref_3x3, ref_brightness_vector);
    max_error = std::max(max_error, *std::max_element(cie_ref3x3_de.cbegin(), cie_ref3x3_de.cend()));
    std::cout << "cie <> ref 3x3: ";
    printDeltaEInfo(cie_ref3x3_de);

    // reference spectra, measured lookup
    image_d image_ref_mL = transformResponse_lookup(ref_response, ref_lookup);
    storeViewableImages(image_ref_mL, ref_brightness_vector, "../spectra_ref_measured_lookup_xyz.ppm", "../spectra_ref_measured_lookup_srgb.ppm", "../spectra_ref_measured_lookup_rec2020.ppm");
    std::vector<double> cie_ref_ml_de = computeImageDeltaE(cie_response, cie_brightness_vector, image_ref_mL, ref_brightness_vector);
    max_error = std::max(max_error, *std::max_element(cie_ref_ml_de.cbegin(), cie_ref_ml_de.cend()));
    std::cout << "cie <> ref mL: ";
    printDeltaEInfo(cie_ref_ml_de);

    // reference spectra, fitted lookup
    image_d image_ref_fL = transformResponse_lookup(ref_response, lookup);
    storeViewableImages(image_ref_fL, ref_brightness_vector, "../spectra_ref_fitted_lookup_xyz.ppm", "../spectra_ref_fitted_lookup_srgb.ppm", "../spectra_ref_fitted_lookup_rec2020.ppm");
    std::vector<double> cie_ref_fL_de = computeImageDeltaE(cie_response, cie_brightness_vector, image_ref_fL, ref_brightness_vector);
    max_error = std::max(max_error, *std::max_element(cie_ref_fL_de.cbegin(), cie_ref_fL_de.cend()));
    std::cout << "cie <> ref fL: ";
    printDeltaEInfo(cie_ref_fL_de);

    // save chromatic error images
    visualizeChromaticError(image_ref_3x3, cie_ref3x3_de, max_error, "../chromatic_cie_ref_3x3.ppm");
    visualizeChromaticError(image_ref_fL, cie_ref_fL_de, max_error, "../chromatic_cie_ref_fL.ppm");
    visualizeChromaticError(image_ref_mL, cie_ref_fL_de, max_error, "../chromatic_cie_ref_mL.ppm");
    visualizeChromaticError(image_fitted_fitted, cie_fitted_fitted_de, max_error, "../chromatic_cie_fitted_fitted.ppm");

    // save error images scaled with max error
    saveErrorImage(image_fitted_fitted.width, image_fitted_fitted.height, max_error, cie_fitted_fitted_de, "../cie__fitted_fitted.ppm");
    saveErrorImage(image_ref_3x3.width, image_ref_3x3.height, max_error, cie_ref3x3_de, "../cie__ref_3x3.ppm");
    saveErrorImage(image_ref_mL.width, image_ref_mL.height, max_error, cie_ref_ml_de, "../cie__ref_measuredL.ppm");
    saveErrorImage(image_ref_fL.width, image_ref_fL.height, max_error, cie_ref_fL_de, "../cie__ref_fittedL.ppm");

    double max_color_err = 0.0000001;
    std::vector<double> cie_ref_3x3_cd = computeImageColorDifference(cie_response, image_ref_3x3);
    max_color_err = std::max(max_color_err, *std::max_element(cie_ref_3x3_cd.cbegin(), cie_ref_3x3_cd.cend()));

    std::vector<double> cie_ref_mL_cd = computeImageColorDifference(cie_response, image_ref_mL);
    max_color_err = std::max(max_color_err, *std::max_element(cie_ref_mL_cd.cbegin(), cie_ref_mL_cd.cend()));

    std::vector<double> cie_ref_fL_cd = computeImageColorDifference(cie_response, image_ref_fL);
    max_color_err = std::max(max_color_err, *std::max_element(cie_ref_fL_cd.cbegin(), cie_ref_fL_cd.cend()));

    std::vector<double> cie_fit_fit_cd = computeImageColorDifference(cie_response, image_fitted_fitted);
    max_color_err = std::max(max_color_err, *std::max_element(cie_fit_fit_cd.cbegin(), cie_fit_fit_cd.cend()));

    std::cout << "ColorDiff: cie <> fit fit: ";
    printColorDeltaInfo(cie_fit_fit_cd);
    std::cout << "ColorDiff: cie <> ref 3x3: ";
    printColorDeltaInfo(cie_ref_3x3_cd);
    std::cout << "ColorDiff: cie <> ref mL: ";
    printColorDeltaInfo(cie_ref_mL_cd);
    std::cout << "ColorDiff: cie <> ref fL: ";
    printColorDeltaInfo(cie_ref_fL_cd);

    saveErrorImage(image_fitted_fitted.width, image_fitted_fitted.height, max_color_err, cie_fit_fit_cd, "../cd_cie__fit_fit.ppm");
    saveErrorImage(image_ref_3x3.width, image_ref_3x3.height, max_color_err, cie_ref_3x3_cd, "../cd_cie__ref_3x3.ppm");
    saveErrorImage(image_ref_mL.width, image_ref_mL.height, max_color_err, cie_ref_mL_cd, "../cd_cie__ref_mL.ppm");
    saveErrorImage(image_ref_fL.width, image_ref_fL.height, max_color_err, cie_ref_fL_cd, "../cd_cie__ref_fL.ppm");

    return 0;
}