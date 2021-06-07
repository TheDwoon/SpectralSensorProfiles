#ifndef SENSOR_SPECTRA_K_MEANS_H
#define SENSOR_SPECTRA_K_MEANS_H

#include <vector>

struct k_means {
    double c1_mean[2];
    std::vector<size_t> c1_idx;
    double c2_mean[2];
    std::vector<size_t> c2_idx;

    double clusterDistance() const;
};

k_means k_means_build(const double* points, size_t count);

#endif //SENSOR_SPECTRA_K_MEANS_H
