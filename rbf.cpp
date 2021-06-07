#include "rbf.h"
#include <iostream>

double rbf::rbf_kernel(double x) const {
    return exp(-m_weight * m_weight * x * x);
}

Eigen::MatrixXd compute_distance_matrix(const Eigen::MatrixX2d &xs) {
    const size_t num_rows = xs.rows();
    Eigen::MatrixXd m(num_rows, num_rows);

    for (Eigen::Index i = 0; i < num_rows; ++i) {
        for (Eigen::Index j = 0; j < num_rows; ++j) {
            double dx = xs(i, 0) - xs(j, 0);
            double dy = xs(i, 1) - xs(j, 1);
            m(i, j) = sqrt(dx * dx + dy * dy);
        }
    }

    return m;
}

rbf::rbf() : m_weight(0) {
}

rbf::rbf(double weight, const double *xs, const double *ys, size_t count) : m_weight(weight) {
    m_xs = Eigen::MatrixX2d();
    m_xs.resize(count, 2);
    for (Eigen::Index i = 0; i < count; ++i) {
        m_xs(i, 0) = xs[i * 2 + 0];
        m_xs(i, 1) = xs[i * 2 + 1];
    }

    Eigen::MatrixX2d fs;
    fs.resize(count, 2);
    for (Eigen::Index i = 0; i < count; ++i) {
        fs(i, 0) = ys[i * 2 + 0];
        fs(i, 1) = ys[i * 2 + 1];
    }

    Eigen::MatrixXd phi = compute_distance_matrix(m_xs);
    for (Eigen::Index i = 0; i < phi.rows(); ++i) {
        for (Eigen::Index j = 0; j < phi.cols(); ++j) {
            phi(i, j) = rbf_kernel(phi(i, j));
        }
    }

    // Validate inverse is correct
    assert(!std::isnan(phi(0, 0)));
    assert(!std::isnan(phi.inverse()(0, 0)));
    assert(!std::isnan(fs(0, 0)));
    m_weights = phi.inverse() * fs;
}

Eigen::VectorXd rbf::eval(const double *x) const {
    Eigen::VectorXd phi;
    phi.resize(m_xs.rows());
    for (Eigen::Index i = 0; i < m_xs.rows(); ++i) {
        double dx = x[0] - m_xs(i, 0);
        double dy = x[1] - m_xs(i, 1);
        double dist = sqrt(dx * dx + dy * dy);
        phi(i) = rbf_kernel(dist);
    }

    return phi.transpose() * m_weights;
}
