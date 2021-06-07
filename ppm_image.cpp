#include "ppm_image.h"
#include <fstream>
#include <iostream>
#include <cassert>

#define PIXEL_IDX(X, Y) (((Y) * m_width + (X)) * 3)

ppm_image::ppm_image() : m_width(0), m_height(0) {}

ppm_image::ppm_image(uint32_t width, uint32_t height): m_width(width), m_height(height) {
    resize(width, height);
}

const uint8_t *ppm_image::read(uint32_t x, uint32_t y) const {
    return m_image.data() + PIXEL_IDX(x, y);
}

void ppm_image::write(uint32_t x, uint32_t y, uint8_t r, uint8_t g, uint8_t b) {
    assert(x >= 0); assert(x < m_width);
    assert(y >= 0); assert(y < m_height);

    uint8_t* ptr = m_image.data() + PIXEL_IDX(x, y);
    ptr[0] = r;
    ptr[1] = g;
    ptr[2] = b;
}

void ppm_image::resize(uint32_t width, uint32_t height) {
    m_width = width;
    m_height = height;
    m_image.resize(3 * m_width * m_height, 0);
}

void ppm_image::load(const std::string &path) {
    std::ifstream file(path, std::ios::in | std::ios::binary);
    assert(file.good());

    char image_type[2] { 0 };
    file.read(image_type, sizeof(image_type));

    // only supported format
    assert(image_type[0] == 'P');
    assert(image_type[1] == '6');

    uint16_t width, height, maxvalue;
    file >> width >> height >> maxvalue;

    char ignored;
    file.read(&ignored, sizeof(char));

    // I don't want to support the bigger file format for now
    assert(maxvalue == 255);

    resize(width, height);

    file.read((char*) m_image.data(), width * height * 3);

    assert(file.good());
    file.close();
}

void ppm_image::store(const std::string &path) {
    std::ofstream file(path, std::ios::out | std::ios::binary | std::ios::trunc);
    assert(file.good());

    file << "P6 " << m_width << " " << m_height << " " << "255\n";
    file.write((const char*)m_image.data(), 3 * m_width * m_height);
    file.flush();

    assert(file.good());
    file.close();
}
