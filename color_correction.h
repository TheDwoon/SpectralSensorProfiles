/*
 * The methods crbt_5f, cbrta_halleyf, lab_f and dt_XYZ_to_lab are taken from
 * https://github.com/darktable-org/darktable/blob/master/src/common/colorspaces_inline_conversions.h#L210
 */

#ifndef SENSOR_SPECTRA_COLOR_CORRECTION_H
#define SENSOR_SPECTRA_COLOR_CORRECTION_H

#include <cmath>

 static inline float cbrt_5f(float f)
{
    uint32_t *p = (uint32_t *)&f;
    *p = *p / 3 + 709921077;
    return f;
}

static inline float cbrta_halleyf(const float a, const float R)
{
    const float a3 = a * a * a;
    const float b = a * (a3 + R + R) / (a3 + a3 + R);
    return b;
}

static inline float lab_f(const float x)
{
    const float epsilon = 216.0f / 24389.0f;
    const float kappa = 24389.0f / 27.0f;
    return (x > epsilon) ? cbrta_halleyf(cbrt_5f(x), x) : (kappa * x + 16.0f) / 116.0f;
}

static inline void dt_XYZ_to_Lab(const float XYZ[3], float Lab[3])
{
    float f[3] = { 0.0f };
    for(int i = 0; i < 3; i++) f[i] = lab_f(XYZ[i]);
    Lab[0] = 116.0f * f[1] - 16.0f;
    Lab[1] = 500.0f * (f[0] - f[1]);
    Lab[2] = 200.0f * (f[1] - f[2]);
}

static inline float compute_delta_e(const float* a, const float* b) {
     float a_lab[3];
     float b_lab[3];

    dt_XYZ_to_Lab(a, a_lab);
    dt_XYZ_to_Lab(b, b_lab);

    float sum { 0.f };
    for (int i = 0; i < 3; ++i)
        sum += (a_lab[i] - b_lab[i]) * (a_lab[i] - b_lab[i]);

    return std::sqrt(sum);
}

#endif //SENSOR_SPECTRA_COLOR_CORRECTION_H

#define GAMMA(X) ((X) <= 0.0031308 ? (X) * 12.92 : std::pow((X), 1.0 / 2.4) * (1 + 0.055) - 0.055)
#define GAMMA_INVERSE(X) ((X) <= 0.04045 ? (X) / 12.92 : std::pow(((X) + 0.055) / (1 + 0.055), 2.4))