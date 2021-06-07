#include "spectral_image.h"
#include <fstream>
#include <cassert>
#include <cmath>
#include <iostream>
#include "color_correction.h"


void spectral_image::load(const std::string &path) {
    std::ifstream file(path, std::ios::in | std::ios::binary);
    assert(file.good());

    // read file header
    file >> m_width >> m_height >> m_dimension;
    char ignore;
    file.read(&ignore, sizeof(char));
    assert(ignore == '\n');

    // read binary content
    m_data.resize(m_width * m_height * m_dimension);
    file.read((char*)m_data.data(), m_width * m_height * m_dimension * sizeof(double));
    assert(file.good());

    file.close();
}

const double *spectral_image::readImage(uint32_t num) const {
    assert(num < m_dimension);

    return m_data.data() + num * m_width * m_height;
}

const double *spectral_image::readPixel(uint32_t dim, uint32_t x, uint32_t y) const {
    const double* image_ptr = readImage(dim);
    return image_ptr + m_width * y + x;
}

void spectral_image::gather(double *buffer, uint32_t x, uint32_t y) const {
    for (uint32_t dim = 0; dim < m_dimension; ++dim) {
        buffer[dim] = *readPixel(dim, x, y);
    }
}

void spectral_image::gammaCorrection() {
    for (size_t i = 0; i < m_data.size(); ++i) {
        if (m_data[i] < 0.)
            std::cout << m_data[i] << std::endl;
        assert(m_data[i] >= 0.);
        m_data[i] = GAMMA_INVERSE(m_data[i]);
    }
}
