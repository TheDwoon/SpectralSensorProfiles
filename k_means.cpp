#include "k_means.h"
#include <cmath>
#include <cassert>
#include <cstdint>

inline double km_distance_squared(const double *p1, const double *p2) {
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    return dx * dx + dy * dy;
}

inline double km_distance(const double *p1, const double *p2) {
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    return sqrt(dx * dx + dy * dy);
}

const double* km_findSmallestPoint(const double* points, size_t count) {
    assert(count >= 1);

    const double* smallest = points;
    for (size_t i = 1; i < count; ++i) {
        if (points[i * 2 + 0] + points[i * 2 + 1] < smallest[0] + smallest[1])
            smallest = points + i * 2;
    }

    return smallest;
}

const double* km_findBiggestPoint(const double* points, size_t count) {
    assert(count >= 1);

    const double* biggest = points;
    for (size_t i = 1; i < count; ++i) {
        if (points[i * 2 + 0] + points[i * 2 + 1] > biggest[0] + biggest[1])
            biggest = points + i * 2;
    }

    return biggest;
}

bool km_cluster_equals(const std::vector<uint8_t> &c1, const std::vector<uint8_t> &c2) {
    assert(c1.size() == c2.size());

    for (size_t i = 0; i < c1.size(); ++i) {
        if (c1[i] != c2[i])
            return false;
    }

    return true;
}

k_means k_means_build(const double *points, size_t count) {
    // seed with min and max element
    // this is probably only good for my use case. But that's what this class is for!
    const double* const smallestPoint = km_findSmallestPoint(points, count);
    const double* const biggestPoint = km_findBiggestPoint(points, count);
    double c1_mean[2] { smallestPoint[0], smallestPoint[1] };
    double c2_mean[2] { biggestPoint[0], biggestPoint[1] };

    std::vector<uint8_t> cluster_idx_old(count);
    std::fill(cluster_idx_old.begin(), cluster_idx_old.end(), 255);

    std::vector<uint8_t> cluster_idx(count);
    std::fill(cluster_idx.begin(), cluster_idx.end(), 254);

    while (!km_cluster_equals(cluster_idx, cluster_idx_old)) {
        // sort points into clusters
        for (size_t i = 0; i < count; ++i) {
            double c1_distance = km_distance_squared(points + 2 * i, c1_mean);
            double c2_distance = km_distance_squared(points + 2 * i, c2_mean);
            if (c1_distance < c2_distance)
                cluster_idx[i] = 1;
            else
                cluster_idx[i] = 2;
        }

        // update cluster mean
        size_t c1_count = 0;
        c1_mean[0] = 0.;
        c1_mean[1] = 0.;

        size_t c2_count = 0;
        c2_mean[0] = 0.;
        c2_mean[1] = 0.;
        for (size_t i = 0; i < count; ++i) {
            if (cluster_idx[i] == 1) {
                c1_mean[0] += points[i * 2 + 0];
                c1_mean[1] += points[i * 2 + 1];
                c1_count += 1;
            } else if (cluster_idx[i] == 2) {
                c2_mean[0] += points[i * 2 + 0];
                c2_mean[1] += points[i * 2 + 1];
                c2_count += 1;
            }
        }

        assert(c1_count > 0);
        assert(c2_count > 0);

        if (c1_count > 0) {
            c1_mean[0] /= c1_count;
            c1_mean[1] /= c1_count;
        }
        if (c2_count > 0) {
            c2_mean[0] /= c2_count;
            c2_mean[1] /= c2_count;
        }

        std::swap(cluster_idx_old, cluster_idx);
    }

    k_means result;
    result.c1_mean[0] = c1_mean[0];
    result.c1_mean[1] = c1_mean[1];
    result.c2_mean[0] = c2_mean[0];
    result.c2_mean[1] = c2_mean[1];

    for (size_t i = 0; i < count; ++i) {
        if (cluster_idx[i] == 1) {
            result.c1_idx.push_back(i);
        } else if (cluster_idx[i] == 2) {
            result.c2_idx.push_back(i);
        }
    }

    return result;
}

double k_means::clusterDistance() const {
    double dx = c1_mean[0] - c2_mean[0];
    double dy = c1_mean[1] - c2_mean[1];
    assert(dx * dx + dy * dy > 0.);
    return sqrt(dx * dx + dy * dy);
}
