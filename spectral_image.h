#ifndef SENSOR_SPECTRA_SPECTRAL_IMAGE_H
#define SENSOR_SPECTRA_SPECTRAL_IMAGE_H


#include <cstdint>
#include <string>
#include <vector>

struct spectral_image {
    uint32_t m_width;
    uint32_t m_height;
    uint32_t m_dimension;
    std::vector<double> m_data;

    void gammaCorrection();
    void gather(double *buffer, uint32_t x, uint32_t y) const;
    void load(const std::string &path);
    inline const double* readImage(uint32_t num) const;
    inline const double* readPixel(uint32_t dim, uint32_t x, uint32_t y) const;
};


#endif //SENSOR_SPECTRA_SPECTRAL_IMAGE_H
