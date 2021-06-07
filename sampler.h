#ifndef SENSOR_SPECTRA_SAMPLER_H
#define SENSOR_SPECTRA_SAMPLER_H

#include <random>
#include <Eigen/Dense>

class Sampler {
private:
    std::default_random_engine randomEngine;
    std::uniform_real_distribution<double> distribution;

public:
    Sampler() {
        std::random_device r;
        distribution = std::uniform_real_distribution<double>(-1, 1);
        randomEngine = std::default_random_engine(r());
    };

    double sample(double base, double delta, int n);
    Eigen::Vector3d sampleVec3(double base, double delta, int n);
};


#endif //SENSOR_SPECTRA_SAMPLER_H
