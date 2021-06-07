#include <iostream>
#include <Eigen/Dense>
#include "camera_xy_lookup.h"
#include "ppm_image.h"
#include "cie1931.h"
#include "color_correction.h"

int main(int argc, char** argv) {
    assert(argc == 2);

    std::string path(argv[1]);
    std::cout << "Loading " << path << std::endl;

    ppm_image input_image;
    input_image.load(path);

    ppm_image output_xyz(input_image.m_width, input_image.m_height);
    ppm_image output_srgb(input_image.m_width, input_image.m_height);
    ppm_image output_rec2020(input_image.m_width, input_image.m_height);

    const Eigen::Vector3d ones {1., 1., 1.};
    const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> xyz2srgb((double*)xyz_to_srgb);
    const Eigen::Vector3d xyz2srg_row_sum = xyz2srgb * ones;
    const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> xyz2rec2020((double*)xyz_to_rec2020);
    const Eigen::Vector3d xyz2rec2020_row_sum = xyz2rec2020 * ones;

    std::vector<double> brightness_vector(input_image.m_height * input_image.m_width);
    for (uint32_t y = 0; y < input_image.m_height; ++y) {
        for (uint32_t x = 0; x < input_image.m_width; ++x) {
            const uint8_t *pixel = input_image.read(x, y);
            brightness_vector[y * input_image.m_width + x] = pixel[0] / 255. + pixel[1] / 255. + pixel[2] / 255.;
        }
    }

    double max_brightness = *std::max_element(brightness_vector.cbegin(), brightness_vector.cend());
    spec::camera_xy_lookup lookup = spec::camera_xy_lookup::loadCameraXYLookUp();
    for (uint32_t y = 0; y < input_image.m_height; ++y) {
        for (uint32_t x = 0; x < input_image.m_width; ++x) {
            const uint8_t *pixel = input_image.read(x, y);

            Eigen::Vector3d rgb { pixel[0] / 255., pixel[1] / 255., pixel[2] / 255. };
            const double brightness = rgb.sum() / 3.;
            rgb /= rgb.sum();

            std::array<double, 2> inverted = lookup.access(rgb(0), rgb(1));
            Eigen::Vector3d xyz {inverted[0], inverted[1], 1. - inverted[0] - inverted[1]};
            if (brightness < 0.035)
                output_xyz.write(x, y, 0, 0, 0);
            else
                output_xyz.write(x, y, 255. * xyz(0), 255. * xyz(1), 255. * xyz(2));

            // 1. / max(row_sum) ensures we do not clip and ruin the color ration due to transformation to srgb
            // will make the color darker though
            Eigen::Vector3d srgb = xyz2srgb * xyz;
            srgb(0) = std::clamp(srgb(0), 0., 1.);
            srgb(1) = std::clamp(srgb(1), 0., 1.);
            srgb(2) = std::clamp(srgb(2), 0., 1.);
            srgb *= (3 * brightness) / max_brightness;
            srgb(0) = GAMMA(srgb(0));
            srgb(1) = GAMMA(srgb(1));
            srgb(2) = GAMMA(srgb(2));
            output_srgb.write(x, y, 255. * srgb(0), 255. * srgb(1), 255. * srgb(2));

            Eigen::Vector3d rec2020 = xyz2rec2020 * xyz * (1. / xyz2rec2020_row_sum.maxCoeff());
            rec2020(0) = std::clamp(rec2020(0), 0., 1.);
            rec2020(1) = std::clamp(rec2020(1), 0., 1.);
            rec2020(2) = std::clamp(rec2020(2), 0., 1.);
            rec2020 /= rec2020.sum();
            rec2020 *= brightness;
            rec2020(0) = GAMMA(rec2020(0));
            rec2020(1) = GAMMA(rec2020(1));
            rec2020(2) = GAMMA(rec2020(2));
            output_rec2020.write(x, y, 255. * rec2020(0), 255. * rec2020(1), 255. * rec2020(2));
        }
    }

    output_xyz.store("../corrected_xyz.ppm");
    output_srgb.store("../corrected_srgb.ppm");
    output_rec2020.store("../corrected_rec2020.ppm");
}
