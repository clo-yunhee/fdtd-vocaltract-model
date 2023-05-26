#ifndef IMAGE_MAKE_WAVE_IMAGE_HH_
#define IMAGE_MAKE_WAVE_IMAGE_HH_

#include "../types.hh"
#include "image.hh"

namespace image {

Image makeWaveImage(const array& frame, double min, double max);

}  // namespace image

#endif  // IMAGE_MAKE_WAVE_IMAGE_HH_