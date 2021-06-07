#ifndef SENSOR_SPECTRA_RBF_H
#define SENSOR_SPECTRA_RBF_H

#include <Eigen/Dense>

struct rbf {
    double m_weight;
    Eigen::MatrixX2d m_xs;
    Eigen::MatrixX2d m_weights;

    rbf();
    rbf(double weight, const double *xs, const double *ys, size_t count);
    Eigen::VectorXd eval(const double* x) const;

private:
    inline double rbf_kernel(double x) const;
};


#endif //SENSOR_SPECTRA_RBF_H
