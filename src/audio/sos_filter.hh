#ifndef AUDIO_SOS_FILTER_HH_
#define AUDIO_SOS_FILTER_HH_

#include <array>
#include <complex>
#include <vector>

namespace audio {

class SOSFilter {
   public:
    SOSFilter(const std::vector<std::array<double, 6>>& sos = {});

    std::vector<double> filter(const std::vector<double>& x);

    const std::vector<std::array<double, 6>>& coefficients() const;

   protected:
    std::vector<std::array<double, 6>> m_sos;
    std::vector<std::array<double, 2>> m_zi;
};

}  // namespace audio

#endif  // AUDIO_SOS_FILTER_HH_