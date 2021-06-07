#include "camera_xy_lookup.h"

#include <iostream>
#include <cassert>
#include <libs/nanoflann.hpp>
#include <clip.h>
#include <numeric>
#include "ppm_image.h"
#include "k_means.h"
#include "rbf.h"
#include "color_correction.h"

spec::camera_xy_lookup::camera_xy_lookup(uint16_t precision) : precision(precision), xStep(1. / (precision - 1)), yStep(1. / (precision - 1)) {
    // ^(\d+)\s([01]\.\d{0,})\s([01]\.\d{0,})\s([01]\.\d{0,})$
}

void spec::camera_xy_lookup::store(const std::string& path) const {
    std::ofstream file(path, std::ios::out | std::ios::trunc | std::ios::binary);
    assert(file.good());

    file.write((char*) &precision, sizeof(precision));
    file.write((char*) raw.data(), sizeof(double) * precision * precision * 2);

    file.flush();
    assert(file.good());

    file.close();
}

spec::camera_xy_lookup spec::camera_xy_lookup::loadCameraXYLookUp(const std::string& path) {
    std::ifstream file(path, std::ios::in | std::ios::binary);
    assert(file.good());

    uint16_t precision { 0 };

    file.read((char*) &precision, sizeof(precision));
    assert(file.good());

    spec::camera_xy_lookup lookup(precision);

    double buffer[2];
    for (uint32_t i = 0; i < precision * precision; ++i) {
        file.read((char*) buffer, sizeof(buffer));

        lookup.raw.push_back(buffer[0]);
        lookup.raw.push_back(buffer[1]);
    }

    assert(file.good());
    file.close();

    lookup.build();
    return lookup;
}

/*
 *  xy -> camera
 *
 *  rg(b) -> xy
 */
spec::camera_xy_lookup
spec::camera_xy_lookup::generateXYLookUp(const spec::Spectrum &red, const spec::Spectrum &green, const spec::Spectrum &blue,
                                         const spec::xy_conversion &table, uint16_t precision) {
    spec::camera_xy_lookup lookup(precision);

    // step through interval [0, 1]
    const double xStep = lookup.xStep;
    const double yStep = lookup.yStep;

    for (uint32_t yIdx = 0; yIdx < precision; ++yIdx) {
        const double y = yIdx * yStep;
        for (uint32_t xIdx = 0; xIdx < precision; ++xIdx) {
            const double x = xIdx * xStep;

            spec::coeff_t coeff = table.access(x, y);
            if (coeff[0] == 0. && coeff[1] == 0. && coeff[2] == 0.) {
                lookup.raw.push_back(-1.);
                lookup.raw.push_back(-1.);
            } else {
                // evaluate input
                spec::Spectrum input(380, 780, coeff.data(), 3);
                double r = spec::interpolatedIntegration(input, red);
                double g = spec::interpolatedIntegration(input, green);
                double b = spec::interpolatedIntegration(input, blue);
                double l = r + g + b;

                lookup.raw.push_back(r / l);
                lookup.raw.push_back(g / l);
            }
        }
    }

    // build lookup
    lookup.filter();
    lookup.build();

    return lookup;
}

struct PointCloud {
    const std::vector<double> &points;

    PointCloud(const std::vector<double> &points) : points(points) {}

    // Must return the number of data poins
    inline size_t kdtree_get_point_count() const { return points.size() / 2; }

    // Must return the dim'th component of the idx'th point in the class:
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return points[2 * idx + dim];
    }

    template<class BBOX>
    inline bool kdtree_get_bbox(BBOX &bb) const {
        return false;
    }
};

// I don't want to type this out more than once :)
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2> my_kd_tree_t;

void spec::camera_xy_lookup::filter() {
    PointCloud cloud(raw);
    my_kd_tree_t index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    ppm_image filtered_elements(precision, precision);

    std::cout << "Filtering..." << std::endl;
    for (uint32_t yIdx = 0; yIdx < precision; ++yIdx) {
        const double y = yIdx * yStep;
        for (uint32_t xIdx = 0; xIdx < precision; ++xIdx) {
            const double x = xIdx * xStep;
            const uint32_t idx = 2 * (yIdx * precision + xIdx);
            const double* rg = &raw[idx];

            if (rg[0] < 0. && rg[1] < 0)
                continue;

            const size_t num_results = 100;
            std::vector<size_t> result_idx(num_results);
            std::vector<double> result_dist(num_results);

            index.knnSearch(rg, num_results, &result_idx[0], &result_dist[0]);
            std::vector<double> points;
            points.reserve(num_results * 2);
            for (size_t i = 0; i < num_results; ++i) {
                points[2 * i + 0] = (result_idx[i] % precision) * xStep;
                points[2 * i + 1] = (result_idx[i] / precision) * yStep;
            }

            k_means cluster = k_means_build(points.data(), num_results);
            if (cluster.clusterDistance() > 0.075) {
                // for some reason it is possible that c1 contains the metamere. Weird but ok!
                const std::vector<size_t> &c_idx =
                        cluster.c1_idx.size() > cluster.c2_idx.size() ? cluster.c2_idx : cluster.c1_idx;
                for (size_t i = 0; i < c_idx.size(); ++i) {
                    const double *p = points.data() + 2 * c_idx[i];
                    uint32_t lx = p[0] / xStep;
                    uint32_t ly = p[1] / yStep;
                    uint32_t wipeIdx = 2 * (ly * precision + lx);

                    filtered_elements.write(raw[wipeIdx + 0] * precision, precision - raw[wipeIdx + 1] * precision,
                                            255, 0, 0);

                    raw[wipeIdx + 0] = -1.;
                    raw[wipeIdx + 1] = -1.;
                }
            }
        }
    }

    filtered_elements.store("../rg_filtered.ppm");
}

#define TARGET_WEIGHT 0.7

void spec::camera_xy_lookup::build() {
    data.resize(2 * precision * precision, -1.);
    assert(raw.size() == 2 * precision * precision);

    PointCloud cloud(raw);
    my_kd_tree_t index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    std::cout << "Building inverse..." << std::endl;

    // Build an RBF for each datapoint we have and interpolate for its neighbors.
    // This will ensure that the matrix will be non singular
    std::vector<rbf> prebuild_rbf(raw.size() / 2);

    const size_t num_interpolate_results = 10;
    for (uint32_t yIdx = 0; yIdx < precision; ++yIdx) {
        for (uint32_t xIdx = 0; xIdx < precision; ++xIdx) {
            const uint32_t global_idx = 2 * (yIdx * precision + xIdx);
            const double point[2] { raw[global_idx + 0], raw[global_idx + 1] };

            if (point[0] < 0. || point[1] < 0.)
                continue;

            // search nearest neighbors used for local interpolation
            std::vector<size_t> result_idx(num_interpolate_results);
            std::vector<double> result_dist(num_interpolate_results);
            index.knnSearch(point, num_interpolate_results, &result_idx[0], &result_dist[0]);

            // collect found points and use them to build an rbf
            std::vector<double> xs;
            xs.reserve(2 * num_interpolate_results);
            std::vector<double> ys;
            ys.reserve(2 * num_interpolate_results);

            for (size_t i = 0; i < num_interpolate_results; ++i) {
                const size_t localIdx = result_idx[i];

                size_t localYIdx = localIdx / precision;
                size_t localXIdx = localIdx % precision;
                double y = localYIdx * yStep;
                double x = localXIdx * xStep;

                if (raw[2 * localIdx + 0] >= 0. && raw[2 * localIdx + 1] >= 0.) {
                    xs.push_back(raw[2 * localIdx + 0]);
                    xs.push_back(raw[2 * localIdx + 1]);
                    ys.push_back(x);
                    ys.push_back(y);
                }
            }

            assert(!xs.empty());
//            double maxDistance { 0 };
//            for (size_t i = 0; i < num_interpolate_results; ++i) {
//                maxDistance += std::sqrt(result_dist[i]);
//            }
//
//            maxDistance /= num_interpolate_results;
//            assert(std::sqrt(*std::max_element(result_dist.begin(), result_dist.end())) < 0.5);
//            double maxDistance = std::sqrt(*std::max_element(result_dist.begin(), result_dist.end()));
            double maxDistance = xStep * 5;
            double weight = std::sqrt(-std::log(TARGET_WEIGHT) / (maxDistance * maxDistance));

            // double maxDistance = std::sqrt(*std::max_element(result_dist.begin(), result_dist.end()));
            prebuild_rbf[yIdx * precision + xIdx] = rbf(weight, xs.data(), ys.data(), xs.size() / 2);
        }
    }

    // actually compute the function inverse using the prebuild rbf
    double* functionInverse = data.data();
    for (uint32_t gIdx = 0; gIdx < precision; ++gIdx) {
        const double g = gIdx * yStep;
        for (uint32_t rIdx = 0; rIdx < precision && rIdx + gIdx < precision; ++rIdx) {
            const double r = rIdx * xStep;
            const uint32_t idx = 2 * (gIdx * precision + rIdx);

            const double point[2] { r, g };

            size_t rbf_index;
            double rbf_distance;

            index.knnSearch(&point[0], 1, &rbf_index, &rbf_distance);

            rbf &rbf = prebuild_rbf[rbf_index];
            Eigen::VectorXd interpolated = rbf.eval(point);

            interpolated(0) = std::clamp(interpolated(0), 0., 1.);
            interpolated(1) = std::clamp(interpolated(1), 0., 1.);
            // make sure the interpolated results are compliant with
            if (interpolated(0) + interpolated(1) > 1.) {
                interpolated /= interpolated.sum();
            }
            assert(!std::isnan(interpolated(0)));
            assert(!std::isnan(interpolated(1)));
            functionInverse[idx + 0] = (rbf_index % precision) * xStep;
            functionInverse[idx + 1] = (rbf_index / precision) * yStep;
//            functionInverse[idx + 0] = interpolated(0);
//            functionInverse[idx + 1] = interpolated(1);
        }
    }
}

std::array<double, 2> spec::camera_xy_lookup::access(double r, double g) const {
    if (r < 0. || g < 0.)
        return {-1., -1.};

    auto xIdx = static_cast<uint32_t>(r / xStep);
    auto yIdx = static_cast<uint32_t>(g / yStep);

    uint32_t gIdx = 2 * (yIdx * precision + xIdx);

    return {data[gIdx + 0] , data[gIdx + 1]};
}

void spec::camera_xy_lookup::storeVectorFile(uint16_t vectorPrecision, const xy_conversion &xy_table,
                                             const spec::Spectrum &cam_red,
                                             const spec::Spectrum &cam_green, const spec::Spectrum &cam_blue,
                                             const std::string &path) const {

    std::ofstream file(path, std::ios::out | std::ios::trunc);
    file << "# x y invertedX invertedY\n";

    std::vector<float> delta_es_img(precision * precision);
    std::vector<float> delta_es;

    for (uint32_t yIdx = 0; yIdx < precision; ++yIdx) {
        const double y = yIdx * xStep;
        for (uint32_t xIdx = 0; xIdx < precision; ++xIdx) {
            const double x = xIdx * yStep;

            // generate spectrum for xyz
            spec::coeff_t cs = xy_table.access(x, y);
            if (cs[0] == 0. && cs[1] == 0. && cs[2] == 0.)
                continue;

            spec::Spectrum spectrum(380, 780, cs.data(), 3);

            Eigen::Vector3d rgb(spec::interpolatedIntegration(spectrum, cam_red),
                                spec::interpolatedIntegration(spectrum, cam_green),
                                spec::interpolatedIntegration(spectrum, cam_blue));
            rgb /= rgb.sum();

            std::array<double, 2> inverted = access(rgb(0), rgb(1));

            if (yIdx % vectorPrecision == 0 && xIdx % vectorPrecision == 0)
                file << x << " " << y << " " << inverted[0] << " " << inverted[1] << "\n";

            const float input_xyz[3] { static_cast<float>(x), static_cast<float>(y), static_cast<float>(1. - x - y) };
            const float output_xyz[3] { static_cast<float>(inverted[0]), static_cast<float>(inverted[1]), static_cast<float>(1. - inverted[0] - inverted[1]) };

            const double de = compute_delta_e(input_xyz, output_xyz);
            delta_es_img[yIdx * precision + xIdx] = de;
            delta_es.push_back(de);
        }
    }

    float min_delta_e = *std::min_element(delta_es.cbegin(), delta_es.cend());
    float max_delta_e = *std::max_element(delta_es.cbegin(), delta_es.cend());
    float avg_delta_e = std::accumulate(delta_es.cbegin(), delta_es.cend(), 0.f) / delta_es.size();

    ppm_image image(precision, precision);
    for (uint32_t yIdx = 0; yIdx < precision; ++yIdx) {
        for (uint32_t xIdx = 0; xIdx < precision; ++xIdx) {
            double de = delta_es_img[yIdx * precision + xIdx];
            de /= max_delta_e;
            image.write(xIdx, precision - 1 - yIdx, 255. * de, 0, 0);
        }
    }
    image.store(path + ".ppm");

    std::sort(delta_es.begin(), delta_es.end());
    size_t clippedSize = 0.9 * delta_es.size();
    float percentile_delta_e = std::accumulate(delta_es.cbegin(), delta_es.cbegin() + clippedSize, 0.f) / clippedSize;
    float median_delta_e = delta_es[delta_es.size() / 2];

    auto it = std::find_if(delta_es.cbegin(), delta_es.cend(), [](float f) { return f > 1.f; });
    size_t distance = std::distance(delta_es.cbegin(), it);

    std::cout << "min \u0394e: " << min_delta_e << " max \u0394e: " << max_delta_e;
    std::cout << " avg \u0394e: " << avg_delta_e << " median \u0394e: " << median_delta_e;
    std::cout << " 90th \u0394e: " << percentile_delta_e << std::endl;
    std::cout << (100 * distance / (float) delta_es_img.size()) << "% elements below 1." << std::endl;
    file.flush();
    file.close();
}

std::array<double, 2> spec::camera_xy_lookup::accessRaw(double x, double y) const {
    if (x < 0. || y < 0. || x + y > 1.)
        return {-1. , -1.};

    auto xIdx = static_cast<uint32_t>(x / xStep);
    auto yIdx = static_cast<uint32_t>(y / yStep);

    uint32_t gIdx = 2 * (yIdx * precision + xIdx);

    return {raw[gIdx + 0] , raw[gIdx + 1]};
}

const double *spec::camera_xy_lookup::readRaw() const {
    return raw.data();
}

