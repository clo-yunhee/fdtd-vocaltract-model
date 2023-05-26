#include "butterworth.hh"

#include <af/constants.h>

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace audio;
using namespace af;

Butterworth::Butterworth() {}

bool Butterworth::loPass(const double fs, const double fc,
                         const uint32_t order) {
    return coefficients(kLoPass, fs, fc, fc, order);
}

bool Butterworth::hiPass(const double fs, const double fc,
                         const uint32_t order) {
    return coefficients(kHiPass, fs, fc, fc, order);
}

bool Butterworth::bandPass(const double fs, const double f1c, const double f2c,
                           const uint32_t order) {
    return coefficients(kBandPass, fs, f1c, f2c, order);
}

bool Butterworth::bandStop(const double fs, const double f1c, const double f2c,
                           const uint32_t order) {
    return coefficients(kBandStop, fs, f1c, f2c, order);
}

bool Butterworth::coefficients(const FilterType filterType, const double fs,
                               const double _f1c, const double _f2c,
                               const uint32_t filterOrder) {
    // ** Pre-warping
    double f1c = 2 * std::tan(Pi * _f1c / fs);
    double f2c = 2 * std::tan(Pi * _f2c / fs);

    // ** Design basic S-plane poles-only analogue LP prototype
    std::vector<dcomplex> poles(2 * filterOrder);

    for (uint32_t k = 0; k < filterOrder; ++k) {
        const double arg = double(2 * k + 1) / double(2 * filterOrder);
        const double real = -std::sin(arg * Pi);
        const double imag = std::cos(arg * Pi);
        poles[2 * k].real(real);
        poles[2 * k].imag(imag);
        poles[2 * k + 1].real(real);
        poles[2 * k + 1].imag(-imag);  // Conjugate.
    }

    std::vector<dcomplex> zeros(0);  // Butterworth LP prototype has no zeros.
    double gain = 1.0;  // Always 1 for the Butterworth prototype lowpass.

    // ** Convert prototype to target filter type (LP/HP/BP/BS) - S-plane

    // Re-orient BP/BS corner frequencies if necessary
    if (f1c > f2c) {
        std::swap(f1c, f2c);
    }

    switch (filterType) {
        case kLoPass:
            convertToLoPass(f2c, zeros, poles, gain);
            break;
        case kHiPass:
            convertToHiPass(f2c, zeros, poles, gain);
            break;
        case kBandPass:
            convertToBandPass(sqrt(f1c * f2c), f2c - f1c, zeros, poles, gain);
            break;
        case kBandStop:
            convertToBandStop(sqrt(f1c * f2c), f2c - f1c, zeros, poles, gain);
            break;
        default:
            std::cerr << "Butterworth: unknown filter type" << std::endl;
            return false;
    }

    // ** SANITY CHECK: Ensure poles are in the left half of the S-plane
    for (const auto& p : poles) {
        if (p.real() > 0) {
            std::cerr << "Butterworth: poles must be in the left half plane"
                      << std::endl;
            return false;
        }
    }

    // ** Map zeros & poles from S-plane to Z-plane
    const double preBltGain = gain;

    for (auto& z : zeros) {
        gain /= bilinear(z);
    }

    for (auto& p : poles) {
        gain *= bilinear(p);
    }

    // ** Split up Z-plane poles & zeros into SOS
    m_sos = zp2sos(zeros, poles);

    // Correct the overall gain
    if (filterType == kLoPass ||
        filterType == kBandPass) {  // pre-BLT is okay for S-plane
        gain = preBltGain *
               (preBltGain / gain);  // 2nd term is how much BLT boosts
    } else if (filterType == kHiPass ||
               filterType == kBandStop) {  // HF gain != DC gain
        gain = 1 / gain;
    }

    for (uint32_t k = 0; k < m_sos.size(); ++k) {
        for (uint32_t i = 0; i < 3; ++i) {
            m_sos[k][i] *= pow(gain, 1 / m_sos.size());
        }
    }

    m_zi.resize(m_sos.size(), {0, 0});

    return true;
}

void Butterworth::convertToLoPass(const double Wc, std::vector<dcomplex>& zeros,
                                  std::vector<dcomplex>& poles, double& gain) {
    gain *= pow(Wc, poles.size());

    zeros.resize(0);         // Poles only
    for (auto& p : poles) {  // Scale poles by the cutoff Wc
        p = Wc * p;
    }
}

void Butterworth::convertToHiPass(const double Wc, std::vector<dcomplex>& zeros,
                                  std::vector<dcomplex>& poles, double& gain) {
    // : hp_S = Wc / lp_S

    // Calculate gain
    dcomplex prodz = 1.0;
    dcomplex prodp = 1.0;

    for (const auto& z : zeros) prodz *= -z;
    for (const auto& p : poles) prodp *= -p;

    gain *= prodz.real() / prodp.real();

    // Convert LP poles to HP
    for (auto& p : poles) {
        if (std::abs(p) != 0) {
            p = Wc / p;
        }
    }

    // Init with zeros, no non-zero values to convert
    zeros.resize(poles.size());
    for (auto& z : zeros) {
        z = 0;
    }
}

void Butterworth::convertToBandPass(const double Wc, const double bw,
                                    std::vector<dcomplex>& zeros,
                                    std::vector<dcomplex>& poles,
                                    double&                gain) {
    // : bp_S = 0.5 * lp_S * BW + 0.5 * sqrt(BW^2 * lp_S^2 - 4*Wc^2)

    // Calculate bandpass gain
    gain *= std::pow(bw, poles.size() - zeros.size());

    const uint32_t numPoles = (uint32_t)poles.size();

    // Convert LP poles to BP in an ordered list of poles and their complex
    // conjugates
    poles.resize(2 * numPoles);

    for (uint32_t i = 0; i < numPoles; ++i) {
        const auto& p = poles[i];
        if (std::abs(p) != 0) {
            const dcomplex term1 = 0.5 * p * bw;
            const dcomplex term2 =
                0.5 * std::sqrt((bw * bw) * (p * p) - (4 * Wc * Wc));
            poles[i] = term1 + term2;
            poles[numPoles + i] = term1 - term2;
        }
    }

    // Init zeros, no non-zero values to convert
    zeros.resize(poles.size());
    for (auto& z : zeros) {
        z = 0;
    }
}

void Butterworth::convertToBandStop(const double Wc, const double bw,
                                    std::vector<dcomplex>& zeros,
                                    std::vector<dcomplex>& poles,
                                    double&                gain) {
    // : bs_S = 0.5 * BW / lp_S + 0.5 * sqrt(BW^2 / lp_S^2 - 4*Wc^2)

    // Calculate gain
    dcomplex prodz = 1.0;
    dcomplex prodp = 1.0;

    for (const auto& z : zeros) prodz *= -z;
    for (const auto& p : poles) prodp *= -p;

    gain *= prodz.real() / prodp.real();

    const uint32_t numPoles = (uint32_t)poles.size();

    // Convert LP zeros to band stop
    zeros.resize(2 * numPoles);

    for (uint32_t i = 0; i < numPoles; ++i) {
        zeros[2 * i].real(0);
        zeros[2 * i].imag(Wc);
        zeros[2 * i + 1].real(0);
        zeros[2 * i + 1].imag(-Wc);  // Complex conjugate
    }

    // Convert LP poles to band stop
    poles.resize(2 * numPoles);

    for (uint32_t i = 0; i < numPoles; ++i) {
        const auto& p = poles[i];
        if (std::abs(p) != 0) {
            const dcomplex term1 = 0.5 * bw / p;
            const dcomplex term2 =
                0.5 * std::sqrt((bw * bw) / (p * p) - (4 * Wc * Wc));
            poles[i] = term1 + term2;
            poles[numPoles + i] = term1 - term2;
        }
    }
}

double Butterworth::bilinear(dcomplex& sz) {
    dcomplex s = sz;
    sz = (2.0 + s) / (2.0 - s);
    // Return the gain
    return std::abs(2.0 - s);
}

std::vector<std::array<double, 6>> Butterworth::zp2sos(
    const std::vector<dcomplex>& zeros, const std::vector<dcomplex>& poles) {
    const uint32_t numZeros = (uint32_t)zeros.size();
    const uint32_t numPoles = (uint32_t)poles.size();

    const uint32_t filterOrder = std::max(numZeros, numPoles);

    // Add zeros at -1, so if S-plane degenerate case
    // where numZeros = 0 will map to -1 in Z-plane.
    std::vector<dcomplex> z(filterOrder, -1.0);
    std::vector<dcomplex> p(filterOrder, 0.0);

    // Copy
    std::copy(zeros.begin(), zeros.end(), z.begin());
    std::copy(poles.begin(), poles.end(), p.begin());

    std::vector<std::array<double, 6>> sos;

    for (uint32_t i = 0; i + 1 < filterOrder; i += 2) {
        sos.push_back({
            1,
            -(z[i] + z[i + 1]).real(),
            (z[i] * z[i + 1]).real(),
            1,
            -(p[i] + p[i + 1]).real(),
            (p[i] * p[i + 1]).real(),
        });
    }

    // Odd filter order thus one pair of poles/zeros remains
    if (filterOrder % 2 == 1) {
        sos.push_back({
            1,
            -z[filterOrder - 1].real(),
            0,
            1,
            -p[filterOrder - 1].real(),
            0,
        });
    }

    return sos;
}