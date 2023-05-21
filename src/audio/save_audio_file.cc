#include "save_audio_file.hh"

#include <sndfile.h>

#include <cstring>

void audio::saveToWavFile(const std::string&         path,
                          const std::vector<double>& audio,
                          const int                  sampleRate) {
    int errnum;

    SF_INFO sfinfo;
    memset(&sfinfo, 0, sizeof(sfinfo));
    sfinfo.samplerate = sampleRate;
    sfinfo.channels = 1;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_DOUBLE;

    SNDFILE* sndfile = sf_open(path.c_str(), SFM_WRITE, &sfinfo);
    if (sndfile == nullptr) {
        fprintf(stderr, "Sndfile open failed: %s\n", sf_strerror(nullptr));
        return;
    }

    sf_count_t framesWritten =
        sf_writef_double(sndfile, audio.data(), audio.size());
    if (audio.size() != framesWritten) {
        fprintf(stderr, "Sndfile didn't write all frames: %ld < %lu\n",
                framesWritten, audio.size());
    }

    errnum = sf_close(sndfile);
    if (errnum != 0) {
        fprintf(stderr, "Sndfile close failed: %s\n", sf_error_number(errnum));
    }
}