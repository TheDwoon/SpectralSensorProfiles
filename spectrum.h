#ifndef SENSOR_SPECTRA_SPECTRUM_H
#define SENSOR_SPECTRA_SPECTRUM_H

#include <glob.h>
#include <cassert>
#include <utility>
#include <vector>
#include <cmath>
#include <string>
#include "rgb2spec.h"

const double redReference[] = {0.0354, 0.0325, 0.0295, 0.0234, 0.0173, 0.0167, 0.0162, 0.0171, 0.0181, 0.016, 0.0139, 0.0132, 0.0125, 0.0116, 0.0107, 0.0108, 0.0109, 0.0107, 0.0105, 0.0116, 0.0126, 0.0122, 0.0118, 0.0156, 0.0195, 0.0296, 0.0397, 0.0655, 0.0912, 0.129, 0.1667, 0.166, 0.1654, 0.157, 0.1487, 0.1652, 0.1817, 0.2369, 0.292, 0.3631, 0.4343, 0.4628, 0.4914, 0.5025, 0.5136, 0.4774, 0.4413, 0.4119, 0.3824, 0.3479, 0.3134, 0.2808, 0.2482, 0.2119, 0.1757, 0.1561, 0.1366, 0.1093, 0.082, 0.0592, 0.0365, 0.023, 0.0095, 0.0086, 0.0076, 0.0073, 0.007, 0.007, 0.007, 0.0061, 0.0052, 0.0042, 0.0033, 0.0024, 0.0015, 0.0017, 0.0019, 0.0021, 0.0023, 0.0025, 0.0027};
const double greenReference[] = {0.0359, 0.0329, 0.03, 0.0248, 0.0197, 0.0218, 0.024, 0.0368, 0.0495, 0.0519, 0.0542, 0.0676, 0.0809, 0.0856, 0.0902, 0.1084, 0.1267, 0.2031, 0.2794, 0.4046, 0.5298, 0.5757, 0.6216, 0.7252, 0.8288, 0.8789, 0.9289, 0.9532, 0.9776, 0.9888, 1., 0.9462, 0.8925, 0.8774, 0.8624, 0.8028, 0.7431, 0.6965, 0.6499, 0.5731, 0.4963, 0.4256, 0.355, 0.2939, 0.2328, 0.1707, 0.1087, 0.0843, 0.0599, 0.0471, 0.0344, 0.0289, 0.0235, 0.019, 0.0146, 0.0137, 0.0128, 0.012, 0.0111, 0.0093, 0.0075, 0.0071, 0.0066, 0.0069, 0.0072, 0.0068, 0.0064, 0.0069, 0.0073, 0.0063, 0.0053, 0.0044, 0.0034, 0.0024, 0.0014, 0.0017, 0.0019, 0.0021, 0.0023, 0.0025, 0.0027};
const double blueReference[] = {0.0334, 0.0327, 0.032, 0.0299, 0.0279, 0.087, 0.1461, 0.3299, 0.5138, 0.5756, 0.6375, 0.6833, 0.7291, 0.7736, 0.818, 0.8347, 0.8515, 0.8425, 0.8336, 0.777, 0.7204, 0.6666, 0.6129, 0.5315, 0.4501, 0.3852, 0.3204, 0.2443, 0.1682, 0.1296, 0.091, 0.0717, 0.0524, 0.0438, 0.0351, 0.0283, 0.0216, 0.0189, 0.0163, 0.0147, 0.0131, 0.0119, 0.0107, 0.0102, 0.0098, 0.0088, 0.0079, 0.0074, 0.0068, 0.0069, 0.0069, 0.0072, 0.0075, 0.0077, 0.0078, 0.008, 0.0082, 0.0076, 0.007, 0.0062, 0.0054, 0.0056, 0.0059, 0.0063, 0.0067, 0.0065, 0.0064, 0.0068, 0.0073, 0.0063, 0.0053, 0.0044, 0.0034, 0.0024, 0.0014, 0.0016, 0.0018, 0.0021, 0.0023, 0.0025, 0.0027};

namespace spec {
    double d_eval_precise(const double* c, const double &lambda, int N);

    class Spectrum {
    public:
        /** start of spectrum_ (included) */
        size_t start;
        /** end of spectrum_ (included) */
        size_t end;
        size_t size;

        std::vector<double> spec;

        Spectrum(size_t start, size_t end, std::vector<double> spec) : start(start), end(end), size(end - start + 1), spec(std::move(spec)) {
            assert(this->spec.size() == size);
        }

        Spectrum(size_t start, size_t end, double* coeff, int n) : start(start), end(end), size(end - start + 1) {
            spec.reserve(size);
            for (size_t i = 0; i < size; ++i) {
                spec.push_back(d_eval_precise(coeff, static_cast<double>(start + i), n));
            }
        }

        Spectrum(size_t start, size_t end, float* coeff) : start(start), end(end), size(end - start + 1) {
            spec.reserve(size);
            for (size_t i = 0; i < size; ++i) {
                spec.push_back(rgb2spec_eval_precise(coeff, static_cast<float>(start + i)));
            }
        }

        Spectrum(size_t start, size_t end, size_t step, const double* data) : start(start), end(end), size(end - start + 1) {
            spec.reserve(size);

            for (size_t i = 0; i < size; ++i) {
                spec.push_back(data[i / step]);
            }
        }

        double getBrightness() const;
        void print(size_t step = 10, float precision = 0.1f);
        void store(const std::string &path) const;
    };

    double integrateSpectrum(const Spectrum &first, const Spectrum &second);
    double interpolatedIntegration(const Spectrum &first, const Spectrum &second);
    Spectrum foldSpectrum(const Spectrum &first, const Spectrum &second);
    void printReferenceSpectra();
    spec::Spectrum fromNormalized(size_t start, size_t end, double* coeff, int n);
    spec::Spectrum loadFromFile(const std::string &path);
    void parseFromFile(spec::Spectrum *red, spec::Spectrum *green, spec::Spectrum *blue, const std::string& path);
}

#endif //SENSOR_SPECTRA_SPECTRUM_H
