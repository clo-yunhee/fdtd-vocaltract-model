#ifndef AUDIO_SAVE_AUDIO_FILE_HH_
#define AUDIO_SAVE_AUDIO_FILE_HH_

#include <string>
#include <vector>

namespace audio {

void saveToWavFile(const std::string& path, const std::vector<double>& audio,
                   const int sampleRate);

}

#endif  // AUDIO_SAVE_AUDIO_FILE_HH_