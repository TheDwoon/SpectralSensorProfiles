#include "spectrum.h"
#include <iostream>
#include <algorithm>
#include <fstream>

void spec::Spectrum::print(size_t step, float precision) {
    const size_t barCount = size / step;
    float* bars = new float[barCount];
    for (size_t lambda = 0; lambda < size; lambda += step) {
        float f = 0;
        for (size_t i = lambda; i < lambda + step; i++)
            f += spec[i];
        f /= step;

        bars[lambda / step] = f / precision;
    }

    int totalHeight = static_cast<int>(1 / precision);
    for (int currentHeight = totalHeight; currentHeight >= 0; currentHeight--) {
        for (size_t i = 0; i < barCount; i++) {
            if (bars[i] >= currentHeight) {
                float left = bars[i] - currentHeight;
                if (left < 1.0 / 8.0)
                    std::cout << "\u2581";
                else if (left < 2.0 / 8.0)
                    std::cout << "\u2582";
                else if (left < 3.0 / 8.0)
                    std::cout << "\u2583";
                else if (left < 4.0 / 8.0)
                    std::cout << "\u2584";
                else if (left < 5.0 / 8.0)
                    std::cout << "\u2585";
                else if (left < 6.0 / 8.0)
                    std::cout << "\u2586";
                else if (left < 7.0 / 8.0)
                    std::cout << "\u2587";
                else
                    std::cout << "\u2588";
            } else
                std::cout << " ";
        }

        std::cout << '\n';
    }
    for (size_t i = 0; i < barCount; i++)
        std::cout << "#";
    std::cout << std::endl;

    delete[] bars;
}

double spec::Spectrum::getBrightness() const {
    return *std::max_element(spec.begin(), spec.end());
}

void spec::Spectrum::store(const std::string &path) const {
    std::ofstream file(path, std::ios::out | std::ios::binary | std::ios::trunc);
    assert(file.good());

    file.write(reinterpret_cast<const char *>(&start), sizeof(start));
    file.write(reinterpret_cast<const char *>(&end), sizeof(end));
    file.write(reinterpret_cast<const char *>(spec.data()), sizeof(double) * spec.size());
    file.flush();
    file.close();
}

double spec::integrateSpectrum(const spec::Spectrum &first, const spec::Spectrum &second) {
    assert(first.size == second.size);
    assert(first.start == second.start);
    assert(first.end == second.end);
    assert(first.spec.size() == first.size);
    assert(first.spec.size() == second.spec.size());

    double acc = 0;
    for (int i = 0; i < first.size; ++i) {
        acc += first.spec[i] * second.spec[i];
    }

    return acc / first.size;
}

double spec::interpolatedIntegration(const spec::Spectrum &first, const spec::Spectrum &second) {
    size_t start = std::max(first.start, second.start);
    size_t end = std::min(first.end, second.end);

    double sum = 0;
    for (size_t lambda = start; lambda <= end; ++lambda) {
        sum += first.spec[lambda - first.start] * second.spec[lambda - second.start];
    }

    return sum / (end - start + 1);
}


spec::Spectrum spec::foldSpectrum(const spec::Spectrum &first, const spec::Spectrum &second) {
    assert(first.size == second.size);
    assert(first.start == second.start);
    assert(first.end == second.end);

    double* data = new double[first.size];
    for (size_t i = 0; i < first.size; ++i) {
        data[i] = first.spec[i] * second.spec[i];
    }

    Spectrum spectrum(first.start, first.end, 1, data);
    delete[] data;

    return spectrum;
}

void spec::printReferenceSpectra() {
    spec::Spectrum redRef(380, 780, 5, (double*) redReference);
    spec::Spectrum greenRef(380, 780, 5, (double *) greenReference);
    spec::Spectrum blueRef(380, 780, 5, (double*) blueReference);

    std::cout << "Reference Red: \n";
    redRef.print(3);
    std::cout << "Reference Green: \n";
    greenRef.print(3);
    std::cout << "Reference Blue: \n";
    blueRef.print(3);
}

double spec::d_eval_precise(const double *c, const double &lambda, int N) {
    // is the same as this just generic and stuff
    //  T x = ((((c1 * lambda + c2) * lambda + c3) * lambda + c4) * lambda + c5) * lambda + c6;
    double x = c[0];
    for (int i = 1; i < N; i++) {
        x = x * lambda + c[i];
    }
    double y = 1. / sqrt(x * x + 1.);
    return 0.5 * x * y + 0.5;
}

spec::Spectrum spec::fromNormalized(size_t start, size_t end, double *coeff, int n) {
    std::vector<double> data;

    size_t size = end - start + 1;
    for (size_t i = 0; i < size; ++i) {
        double lambda = static_cast<double>(i) / (size - 1);
        data.push_back(d_eval_precise(coeff, lambda, n));
    }

    return spec::Spectrum(start, end, data);
}

spec::Spectrum spec::loadFromFile(const std::string &path) {
    std::ifstream file(path, std::ios::in | std::ios::binary);
    assert(file.good());

    size_t start, end;
    file.read(reinterpret_cast<char *>(&start), sizeof(start));
    file.read(reinterpret_cast<char *>(&end), sizeof(end));
    double* data = new double[end - start];
    file.read(reinterpret_cast<char *>(data), sizeof(double) * (end - start + 1));

    spec::Spectrum spectrum(start, end, 1, data);
    delete[] data;

    return spectrum;
}

void spec::parseFromFile(spec::Spectrum *red, spec::Spectrum *green, spec::Spectrum *blue, const std::string& path) {
    std::ifstream file(path, std::ios::in);
    assert(file.good());

    size_t lambda;
    double r, g, b;

    std::vector<size_t> lambdas;
    std::vector<double> r_values;
    std::vector<double> g_values;
    std::vector<double> b_values;

    while (file >> lambda >> r >> g >> b) {
        lambdas.push_back(lambda);
        r_values.push_back(r);
        g_values.push_back(g);
        b_values.push_back(b);
    }

    const size_t min_lambda = lambdas[0];
    const size_t max_lambda = lambdas[lambdas.size() - 1];
    const size_t lambda_step = lambdas[1] - lambdas[0];

    *red = spec::Spectrum(min_lambda, max_lambda + lambda_step - 1, lambda_step, r_values.data());
    *green = spec::Spectrum(min_lambda, max_lambda + lambda_step - 1, lambda_step, g_values.data());
    *blue = spec::Spectrum(min_lambda, max_lambda + lambda_step - 1, lambda_step, b_values.data());
}
