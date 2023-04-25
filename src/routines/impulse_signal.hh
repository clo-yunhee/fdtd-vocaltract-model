#ifndef ROUTINES_IMPULSE_SIGNAL_HH_
#define ROUTINES_IMPULSE_SIGNAL_HH_

#include "../eigen.hh"

namespace src {

// sampleRate    : audio sampling rate
// lcf           : low cut-off frequency
// hcf           : high cut-off frequency
// maxExcitation : excitation magnitude
// [returns]     : excitation velocity
ArrayXd ImpulseSignal(uint32_t sampleRate, double maxExcitation, double lcf,
                      double hcf);

}  // namespace src

#endif  // ROUTINES_IMPULSE_SIGNAL_HH_