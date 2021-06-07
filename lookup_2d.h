#ifndef SENSOR_SPECTRA_LOOKUP_2D_H
#define SENSOR_SPECTRA_LOOKUP_2D_H

#include <vector>
#include <cstdint>
#include <iostream>

class lookup_2d {
private:
    uint32_t m_precision;
    std::vector<double> m_data;

public:
    explicit lookup_2d(uint32_t precision = 512);
    const double* read(double x, double y) const;
    double* write(double x, double y);

    friend std::ostream& operator<<(std::ostream& stream, const lookup_2d& lookup) {
        stream.write((char*) &lookup.m_precision, sizeof(lookup.m_precision));
        stream.write((char*) lookup.m_data.data(), sizeof(double) * lookup.m_precision * lookup.m_precision * 2);
        return stream;
    }
};

#endif //SENSOR_SPECTRA_LOOKUP_2D_H
