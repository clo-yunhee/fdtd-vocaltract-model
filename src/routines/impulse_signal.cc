#include "impulse_signal.hh"

#include "../constants.hh"

using namespace vt::constants;

ArrayXd src::ImpulseSignal(const uint32_t sampleRate,
                           const double maxExcitation, double lcf, double hcf) {
    auto excitationX = zeros<ArrayXd>(sampleRate);

    // Define your cut-off frequency as a fraction of sampling rate
    lcf = lcf / sampleRate;
    hcf = hcf / sampleRate;

    // roll-off: Steepness of transmission frequency between the pass band and
    // stop band
    const double   bw = 0.001;
    const uint32_t M = 4 / bw;  // Number of points

    auto imp = zeros<ArrayXd>(M + 1);
    auto lp_sinc = zeros<ArrayXd>(M + 1);
    auto hp_sinc = zeros<ArrayXd>(M + 1);

    // STEP1: Design a low pass filter using sinc function
    double sumVal = 0;
    for (uint32_t i = 1; i <= M + 1; ++i) {
        if (i - 1 != M / 2) {
            lp_sinc(i) =
                sin(2 * pi * hcf * ((i - 1) - M / 2)) / ((i - 1) - M / 2);
        } else {
            lp_sinc(i) = 2 * pi * hcf;
        }
        // Multiply by Hamming window
        lp_sinc(i) *= (0.54 - 0.46 * cos(2 * pi * (i - 1) / M));
        sumVal += lp_sinc(i);
    }

    // STEP2: Normalize the value of the low pass filter
    for (uint32_t i = 1; i <= M + 1; ++i) {
        lp_sinc(i) /= sumVal;
    }

    // STEP3: Design a high pass filter using sinc function
    // [Note]: We can design a high pass filter using low pass filter and
    // spectral inversion.
    // First designed the low pass filter using cut-off frequency.
    // Then invert the low pass filter and add 1 to the sample at the shape
    // of symmetry. (Chapter 14: The scientist and engineer's guide to signal
    // processing)
    sumVal = 0;
    for (uint32_t i = 1; i <= M + 1; ++i) {
        if (i - 1 != M / 2) {
            hp_sinc(i) =
                sin(2 * pi * lcf * ((i - 1) - M / 2)) / ((i - 1) - M / 2);
        } else {
            hp_sinc(i) = 2 * pi * lcf;
        }
        // Multiply by Hamming window
        hp_sinc(i) *= (0.54 - 0.46 * cos(2 * pi * (i - 1) / M));
        sumVal += hp_sinc(i);
    }

    // Normalize and invert, to turn into high pass sinc
    for (uint32_t i = 1; i <= M + 1; ++i) {
        hp_sinc(i) /= sumVal;
        hp_sinc(i) = -hp_sinc(i);
    }
    hp_sinc(M / 2 + 1) += 1;

    // STEP4: Construct a window sync (bandpass) filter
    for (uint32_t i = 1; i <= M + 1; ++i) {
        imp(i) = hp_sinc(i) + lp_sinc(i);
        imp(i) = -imp(i);
    }
    imp(M / 2 + 1) += 1;

    const auto ex = imp * maxExcitation;

    excitationX(seq(1, M + 1)) = ex(seq(1, M + 1));

    return excitationX;
}
