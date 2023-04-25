#include "generate_audio.hh"

#include "butterworth.hh"

std::vector<double> audio::generateFromPressure(
    const Ref<const ArrayXd> Pr_Audio, const double srateBase,
    const double srate, const uint32_t srateMul) {
    const double maxPressure = Pr_Audio.maxCoeff();

    ArrayXd audioOutput = Pr_Audio / maxPressure;
    audioOutput = (audioOutput > +1).select(+1, audioOutput);
    audioOutput = (audioOutput < -1).select(-1, audioOutput);

    // Generate a low pass butterworth filter (IIR)
    const double fc = srateBase / 2;  // Cutoff frequency
    Butterworth  butter;
    butter.loPass(srate, fc, 2);

    // Filter the original signal
    std::vector<double> filteredAudio(&audioOutput(1), &audioOutput(last));
    filteredAudio = butter.filter(filteredAudio);

    // Downsample the signal
    std::vector<double> finalAudioOutput;
    for (uint32_t i = 0; i < filteredAudio.size(); i += srateMul) {
        finalAudioOutput.push_back(filteredAudio[i]);
    }

    return finalAudioOutput;
}
