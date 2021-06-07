#ifndef SENSOR_SPECTRA_XY_CONVERSION_H
#define SENSOR_SPECTRA_XY_CONVERSION_H

#include <cstdint>
#include <fstream>
#include <vector>
#include <array>

namespace spec {
    typedef std::array<double, 3> coeff_t;

    struct xy_conversion {
        uint32_t magic;
        uint16_t version;
        uint8_t  channels;
        uint8_t  datatype;
        uint32_t width;
        uint32_t height;
        std::vector<float> data;

        [[nodiscard]] coeff_t access(double x, double y) const;
    };

    xy_conversion loadConversion(const std::string& path = "resources/spectra.lut");
}

#endif //SENSOR_SPECTRA_XY_CONVERSION_H
