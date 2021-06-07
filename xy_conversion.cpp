#include "xy_conversion.h"
#include <cassert>
#include <iostream>

#define CONVERSION_HEADER_SIZE 16

spec::xy_conversion spec::loadConversion(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    assert(!file.fail());

    xy_conversion conversion;
    file.read(reinterpret_cast<char *>(&conversion), CONVERSION_HEADER_SIZE);
    assert(file.good());

    assert(conversion.channels == 4);
    assert(conversion.version == 2);

    float buffer[4];
    for (size_t i = 0; i < conversion.width * conversion.height; ++i) {
        file.read(reinterpret_cast<char *>(buffer), sizeof(buffer));
        assert(file.good());

        conversion.data.push_back(buffer[0]);
        conversion.data.push_back(buffer[1]);
        conversion.data.push_back(buffer[2]);
        conversion.data.push_back(buffer[3]);
    }

    file.close();

    return std::move(conversion);
}

#define CLAMP(X, MIN, MAX) ((X) < (MIN) ? (MIN) : ((X) > (MAX) ? (MAX) : (X)))
spec::coeff_t spec::xy_conversion::access(double x, double y) const {
    int xi = CLAMP(x * width + 0.5f, 0, width - 1);
    int yi = CLAMP(y * height + 0.5f, 0, height - 1);

    int idx = 4 * (yi * width + xi);
    assert(idx >= 0 && idx + 2 < data.size());

    spec::coeff_t res = { data[idx], data[idx + 1], data[idx + 2] };
    return res;
}
