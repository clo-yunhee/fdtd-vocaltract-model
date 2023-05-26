#ifndef AUDIO_GENERATE_AUDIO_HH_
#define AUDIO_GENERATE_AUDIO_HH_

#include <vector>

#include "../types.hh"

namespace audio {

std::vector<double> generateFromPressure(const array& Pr_Audio,
                                         double srateBase, double srate,
                                         uint32_t srateMul);

}

#endif  // AUDIO_GENERATE_AUDIO_HH_