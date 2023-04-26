#include "generate_audio.hh"

#include "butterworth.hh"

std::vector<double> audio::generateFromPressure(
    const Ref<const ArrayXd> Pr_Audio, const double srateBase,
    const double srate, const uint32_t srateMul) {
    std::vector<double> audioOutput(length(Pr_Audio));

    double maxPressure = -DBL_MAX;

    for (uint32_t i = 0; i < audioOutput.size(); ++i) {
        audioOutput[i] = Pr_Audio(i + 1);
        if (std::abs(audioOutput[i]) > maxPressure) {
            maxPressure = std::abs(audioOutput[i]);
        }
    }

    for (uint32_t i = 0; i < audioOutput.size(); ++i) {
        audioOutput[i] /= maxPressure;
        if (audioOutput[i] < -1) audioOutput[i] = -1;
        if (audioOutput[i] > +1) audioOutput[i] = +1;
    }

    return audioOutput;

    // Generate a low pass butterworth filter (IIR)
    const double fc = srateBase / 2;  // Cutoff frequency
    Butterworth  butter;
    butter.loPass(srate, fc, 2);

    // Filter the original signal
    std::vector<double> filteredAudio = butter.filter(audioOutput);

    // Downsample the signal
    std::vector<double> finalAudioOutput;
    for (uint32_t i = 0; i < filteredAudio.size(); i += srateMul) {
        finalAudioOutput.push_back(filteredAudio[i]);
    }

    return finalAudioOutput;
}
