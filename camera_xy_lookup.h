#ifndef SENSOR_SPECTRA_CAMERA_XY_LOOKUP_H
#define SENSOR_SPECTRA_CAMERA_XY_LOOKUP_H

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include "spectrum.h"
#include "xy_conversion.h"

namespace spec {
    class camera_xy_lookup {
    private:
        uint16_t precision;
        /** Contains raw function data [xy] -> [rg(b)] */
        std::vector<double> raw;
        /** Contains the inverted function represented by the raw data points stored in raw.
         * [rg(b)] -> [xy]
         */
        std::vector<double> data;

        const double xStep;
        const double yStep;

        /**
         * Filters raw data and removes metameres.
         */
        void filter();
        /**
         * Builds the inverse function from the given raw data.
         */
        void build();
    public:
        explicit camera_xy_lookup(uint16_t precision);

        const double* readRaw() const;

        /**
         * Access raw function data.
         * [xy] -> [rg(b)] with b = 1 - r - g.
         *
         * @param x
         * @param y
         * @return
         */
        std::array<double, 2> accessRaw(double x, double y) const;

        /**
         * Access inverted function at r, g values.
         * [rg] -> [xy]
         *
         * @param r red value of camera
         * @param g green value of camera
         * @return xy used to generate the spectrum resulting in the specific rg values in camera RGB
         */
        std::array<double, 2> access(double r, double g) const;

        void store(const std::string& path = "../camera_lookup.bin") const;
        void storeVectorFile(uint16_t vectorPrecision, const xy_conversion &xy_table, const spec::Spectrum &cam_red,
                             const spec::Spectrum &cam_green, const spec::Spectrum &cam_blue,
                             const std::string &path = "../inverse.dat") const;

        static camera_xy_lookup loadCameraXYLookUp(const std::string& path = "../camera_lookup.bin");
        static camera_xy_lookup generateXYLookUp(const spec::Spectrum &red, const spec::Spectrum &green, const spec::Spectrum &blue,
                                          const xy_conversion &table, uint16_t precision = 1024);
    };
}

#endif //SENSOR_SPECTRA_CAMERA_XY_LOOKUP_H
