#include "impulse_signal.hh"

#include "../constants.hh"

using namespace vt::constants;

array src::ImpulseSignal(const uint32_t sampleRate, const double maxExcitation,
                         double lcf, double hcf) {
    array excitationX = constant(0, sampleRate, 1);

    // Define your cut-off frequency as a fraction of sampling rate
    lcf = lcf / sampleRate;
    hcf = hcf / sampleRate;

    // roll-off: Steepness of transmission frequency between the pass band and
    // stop band
    const double   bw = 0.001;
    const uint32_t M = 4 / bw;  // Number of points

    array imp = constant(0, M + 1, 1);
    array lp_sinc = constant(0, M + 1, 1);
    array hp_sinc = constant(0, M + 1, 1);

    // STEP1: Design a low pass filter using sinc function
    for (uint32_t i = 0; i < M + 1; ++i) {
        if (i != M / 2) {
            lp_sinc(i) = sin(2 * Pi * hcf * (i - .5 * M)) / (i - .5 * M);
        } else {
            lp_sinc(i) = 2 * Pi * hcf;
        }
        // Multiply by Hamming window
        lp_sinc(i) *= (0.54 - 0.46 * cos(2 * Pi * i / M));
    }

    // STEP2: Normalize the value of the low pass filter
    lp_sinc /= sum<float>(lp_sinc);

    // STEP3: Design a high pass filter using sinc function
    // [Note]: We can design a high pass filter using low pass filter and
    // spectral inversion.
    // First designed the low pass filter using cut-off frequency.
    // Then invert the low pass filter and add 1 to the sample at the shape
    // of symmetry. (Chapter 14: The scientist and engineer's guide to signal
    // processing)
    for (uint32_t i = 0; i < M + 1; ++i) {
        if (i != M / 2) {
            hp_sinc(i) = sin(2 * Pi * lcf * (i - .5 * M)) / (i - .5 * M);
        } else {
            hp_sinc(i) = 2 * Pi * lcf;
        }
        // Multiply by Hamming window
        hp_sinc(i) *= (0.54 - 0.46 * cos(2 * Pi * i / M));
    }

    // Normalize and invert, to turn into high pass sinc
    hp_sinc /= sum<float>(hp_sinc);
    hp_sinc = -hp_sinc;
    hp_sinc(M / 2) = hp_sinc(M / 2) + 1;

    // STEP4: Construct a window sync (bandpass) filter
    imp = hp_sinc + lp_sinc;
    imp = -imp;
    imp(M / 2) = imp(M / 2) + 1;

    array ex = imp * maxExcitation;

    if (M + 1 >= sampleRate) {
        excitationX = ex;
    } else {
        excitationX(seq(0, M)) = ex;
    }

    return excitationX;
}
