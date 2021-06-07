#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include "ppm_image.h"

int main(int argc, char** argv) {
    assert(argc == 2);

    std::string path(argv[1]);
    std::cout << "Generating histogram for " << path << std::endl;

    ppm_image image;
    image.load(path);

    std::ofstream file("../image_spectra.dat", std::ios::out | std::ios::trunc);
    assert(file.good());

    std::unordered_set<size_t> hashes;
    std::hash<double> hf;

    for (uint32_t y = 0; y < image.m_height; ++y) {
        for (uint32_t x = 0; x < image.m_width; ++x) {
            const uint8_t* pixel = image.read(x, y);
            double r = pixel[0] / 255.;
            double g = pixel[1] / 255.;
            double b = pixel[2] / 255.;
            double l = r + g + b;

            r /= l;
            g /= l;

            size_t r_hash = hf(r);
            size_t g_hash = hf(g);
            size_t hash = r_hash ^ g_hash;
            if (hashes.find(hash) != hashes.end() || l < 0.05)
                continue;

            hashes.insert(hash);
            file << r << " " << g << '\n';
        }
    }

    file.flush();
    file.close();
}