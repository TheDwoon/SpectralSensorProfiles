/**
 * This file contains:
 *
 * 1. CIE 1931 curves at sampled at 5nm intervals
 *
 * 2. CIE D65 and D50 spectra sampled at 5nm intervals.
 *    Both are normalized to have unit luminance.
 *
 * 3. XYZ <-> sRGB conversion matrices
 *    XYZ <-> ProPhoto RGB conversion matrices
 *
 * 4. A convenience function "cie_interp" to access the discretized
 *    data at arbitrary wavelengths (with linear interpolation)
}
 */
#pragma once
#include "cie1931_c.h"
#define N(x) (x / 10566.864005283874576)

const double cie_d65[CIE_SAMPLES] = {
    N(46.6383),  N(49.3637),  N(52.0891),  N(51.0323),  N(49.9755),  N(52.3118),  N(54.6482),  N(68.7015),
    N(82.7549),  N(87.1204),  N(91.486),   N(92.4589),  N(93.4318),  N(90.057),   N(86.6823),  N(95.7736),
    N(104.865),  N(110.936),  N(117.008),  N(117.41),   N(117.812),  N(116.336),  N(114.861),  N(115.392),
    N(115.923),  N(112.367),  N(108.811),  N(109.082),  N(109.354),  N(108.578),  N(107.802),  N(106.296),
    N(104.79),   N(106.239),  N(107.689),  N(106.047),  N(104.405),  N(104.225),  N(104.046),  N(102.023),
    N(100.0),    N(98.1671),  N(96.3342),  N(96.0611),  N(95.788),   N(92.2368),  N(88.6856),  N(89.3459),
    N(90.0062),  N(89.8026),  N(89.5991),  N(88.6489),  N(87.6987),  N(85.4936),  N(83.2886),  N(83.4939),
    N(83.6992),  N(81.863),   N(80.0268),  N(80.1207),  N(80.2146),  N(81.2462),  N(82.2778),  N(80.281),
    N(78.2842),  N(74.0027),  N(69.7213),  N(70.6652),  N(71.6091),  N(72.979),   N(74.349),   N(67.9765),
    N(61.604),   N(65.7448),  N(69.8856),  N(72.4863),  N(75.087),   N(69.3398),  N(63.5927),  N(55.0054),
    N(46.4182),  N(56.6118),  N(66.8054),  N(65.0941),  N(63.3828),  N(63.8434),  N(64.304),   N(61.8779),
    N(59.4519),  N(55.7054),  N(51.959),   N(54.6998),  N(57.4406),  N(58.8765),  N(60.3125)
};

#undef N

#define N(x) (x / 106.8)
const double cie_e[CIE_SAMPLES] = {
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0),  N(1.0),
    N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0), N(1.0)
};
#undef N

#define N(x) (x / 10503.2)

const double cie_d50[CIE_SAMPLES] = {
    N(23.942000),  N(25.451000),  N(26.961000),  N(25.724000),  N(24.488000),
    N(27.179000),  N(29.871000),  N(39.589000),  N(49.308000),  N(52.910000),
    N(56.513000),  N(58.273000),  N(60.034000),  N(58.926000),  N(57.818000),
    N(66.321000),  N(74.825000),  N(81.036000),  N(87.247000),  N(88.930000),
    N(90.612000),  N(90.990000),  N(91.368000),  N(93.238000),  N(95.109000),
    N(93.536000),  N(91.963000),  N(93.843000),  N(95.724000),  N(96.169000),
    N(96.613000),  N(96.871000),  N(97.129000),  N(99.614000),  N(102.099000),
    N(101.427000), N(100.755000), N(101.536000), N(102.317000), N(101.159000),
    N(100.000000), N(98.868000),  N(97.735000),  N(98.327000),  N(98.918000),
    N(96.208000),  N(93.499000),  N(95.593000),  N(97.688000),  N(98.478000),
    N(99.269000),  N(99.155000),  N(99.042000),  N(97.382000),  N(95.722000),
    N(97.290000),  N(98.857000),  N(97.262000),  N(95.667000),  N(96.929000),
    N(98.190000),  N(100.597000), N(103.003000), N(101.068000), N(99.133000),
    N(93.257000),  N(87.381000),  N(89.492000),  N(91.604000),  N(92.246000),
    N(92.889000),  N(84.872000),  N(76.854000),  N(81.683000),  N(86.511000),
    N(89.546000),  N(92.580000),  N(85.405000),  N(78.230000),  N(67.961000),
    N(57.692000),  N(70.307000),  N(82.923000),  N(80.599000),  N(78.274000),
    N(0),          N(0),          N(0),          N(0),          N(0),
    N(0),          N(0),          N(0),          N(0)
};

#undef N

#define N(x) (x / 10536.3)

const double cie_d60[CIE_SAMPLES] = {
    N(38.683115),  N(41.014457),  N(42.717548),  N(42.264182),  N(41.454941),
    N(41.763698),  N(46.605319),  N(59.226938),  N(72.278594),  N(78.231500),
    N(80.440600),  N(82.739580),  N(82.915027),  N(79.009168),  N(77.676264),
    N(85.163609),  N(95.681274),  N(103.267764), N(107.954821), N(109.777964),
    N(109.559187), N(108.418402), N(107.758141), N(109.071548), N(109.671404),
    N(106.734741), N(103.707873), N(103.981942), N(105.232199), N(105.235867),
    N(104.427667), N(103.052881), N(102.522934), N(104.371416), N(106.052671),
    N(104.948900), N(103.315154), N(103.416286), N(103.538599), N(102.099304),
    N(100.000000), N(97.992725),  N(96.751421),  N(97.102402),  N(96.712823),
    N(93.174457),  N(89.921479),  N(90.351933),  N(91.999793),  N(92.384009),
    N(92.098710),  N(91.722859),  N(90.646003),  N(88.327552),  N(86.526483),
    N(87.034239),  N(87.579186),  N(85.884584),  N(83.976140),  N(83.743140),
    N(84.724074),  N(86.450818),  N(87.493491),  N(86.546330),  N(83.483070),
    N(78.268785),  N(74.172451),  N(74.275184),  N(76.620385),  N(79.423856),
    N(79.051849),  N(71.763360),  N(65.471371),  N(67.984085),  N(74.106079),
    N(78.556612),  N(79.527120),  N(75.584935),  N(67.307163),  N(55.275106),
    N(49.273538),  N(59.008629),  N(70.892412),  N(70.950115),  N(67.163996),
    N(67.445480),  N(68.171371),  N(66.466636),  N(62.989809),  N(58.067786),
    N(54.990892),  N(56.915942),  N(60.825601),  N(62.987850)
};

#undef N

const double xyz_to_srgb[3][3] = {
    { 3.240479, -1.537150, -0.498535 },
    {-0.969256,  1.875991,  0.041556 },
    { 0.055648, -0.204043,  1.057311 }
};

const double srgb_to_xyz[3][3] = {
    { 0.412453, 0.357580, 0.180423 },
    { 0.212671, 0.715160, 0.072169 },
    { 0.019334, 0.119193, 0.950227 }
};

const double xyz_to_xyz[3][3] = {
  { 1.0, 0.0, 0.0 },
  { 0.0, 1.0, 0.0 },
  { 0.0, 0.0, 1.0 },
};

const double xyz_to_ergb[3][3] = {
  {  2.689989, -1.276020, -0.413844},
  { -1.022095,  1.978261,  0.043821},
  {  0.061203, -0.224411,  1.162859},
};

const double ergb_to_xyz[3][3] = {
  { 0.496859,  0.339094,  0.164047 },
  { 0.256193,  0.678188,  0.065619 },
  { 0.023290,  0.113031,  0.863978 },
};

const double xyz_to_prophoto_rgb[3][3] = {
    { 1.3459433,  -0.2556075, -0.0511118 },
    { -0.5445989,  1.5081673,  0.0205351 },
    {  0.0000000,  0.0000000,  1.2118128 }
};

const double prophoto_rgb_to_xyz[3][3] = {
    { 0.7976749,  0.1351917,  0.0313534 },
    { 0.2880402,  0.7118741,  0.0000857 },
    { 0.0000000,  0.0000000,  0.8252100 }
};

const double xyz_to_aces2065_1[3][3] = {
    {  1.0498110175, 0.0000000000, -0.0000974845 },
    { -0.4959030231, 1.3733130458, 0.0982400361 },
    {  0.0000000000, 0.0000000000, 0.9912520182 }
};

const double aces2065_1_to_xyz[3][3] = {
    { 0.9525523959, 0.0000000000, 0.0000936786 },
    { 0.3439664498, 0.7281660966, -0.0721325464 },
    { 0.0000000000, 0.0000000000, 1.0088251844 }
};

const double xyz_to_rec2020[3][3] = {
    {  1.7166511880, -0.3556707838, -0.2533662814 },
    { -0.6666843518,  1.6164812366,  0.0157685458 },
    {  0.0176398574, -0.0427706133,  0.9421031212 }
};

const double rec2020_to_xyz[3][3] = {
    { 0.6369580483, 0.1446169036, 0.1688809752 },
    { 0.2627002120, 0.6779980715, 0.0593017165 },
    { 0.0000000000, 0.0280726930, 1.0609850577 }
};

double cie_interp(const double *data, double x) {
    x -= CIE_LAMBDA_MIN;
    x *= (CIE_SAMPLES - 1) / (CIE_LAMBDA_MAX - CIE_LAMBDA_MIN);
    int offset = std::min(std::max(0, (int) x), CIE_SAMPLES - 2);
    double weight = x - offset;
    return (1.0 - weight) * data[offset] + weight * data[offset + 1];
}
