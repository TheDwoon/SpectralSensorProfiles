#ifndef SENSOR_SPECTRA_PPM_IMAGE_H
#define SENSOR_SPECTRA_PPM_IMAGE_H

#include <vector>
#include <cstdint>
#include <string>

struct ppm_image {
    uint32_t m_width;
    uint32_t m_height;
    std::vector<uint8_t> m_image;

    ppm_image();
    ppm_image(uint32_t width, uint32_t height);

    const uint8_t* read(uint32_t x, uint32_t y) const;
    void write(uint32_t x, uint32_t y, uint8_t r, uint8_t g, uint8_t b);
    void resize(uint32_t width, uint32_t height);

    void load(const std::string &path);
    void store(const std::string &path);
};


#endif //SENSOR_SPECTRA_PPM_IMAGE_H
