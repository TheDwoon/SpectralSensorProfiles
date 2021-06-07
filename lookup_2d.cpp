#include "lookup_2d.h"

lookup_2d::lookup_2d(uint32_t precision) : m_precision(precision), m_data(std::vector<double>(2 * precision)) {

}

const double *lookup_2d::read(double x, double y) const {
    if (x < 0. || x > 1. || y < 0. || y > 1.)
        return nullptr;

    const double step = 1. / (m_precision - 1);
    const uint32_t xIdx = static_cast<uint32_t>(x / step);
    const uint32_t yIdx = static_cast<uint32_t>(y / step);
    const uint32_t idx = 2 * (yIdx * m_precision + xIdx);

    return m_data.data() + idx;
}

double *lookup_2d::write(double x, double y) {
    if (x < 0. || x > 1. || y < 0. || y > 1.)
        return nullptr;

    const double step = 1. / (m_precision - 1);
    const uint32_t xIdx = static_cast<uint32_t>(x / step);
    const uint32_t yIdx = static_cast<uint32_t>(y / step);
    const uint32_t idx = 2 * (yIdx * m_precision + xIdx);

    return m_data.data() + idx;
}
