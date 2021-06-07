#include "sampler.h"

double Sampler::sample(double base, double delta, int n) {
    double s = base;
    for (int i = 0; i < n; ++i) {
        s += distribution(randomEngine) * delta;
    }

    return s;
}

Eigen::Vector3d Sampler::sampleVec3(double base, double delta, int n) {
    double x = sample(base, delta, n);
    double y = sample(base, delta, n);
    double z = sample(base, delta, n);

    return {x, y, z};
}
