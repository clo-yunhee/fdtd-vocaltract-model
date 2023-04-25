#ifndef IMAGE_SAVE_TO_PNG_HH_
#define IMAGE_SAVE_TO_PNG_HH_

#include <string>

#include "image.hh"

namespace image {

void saveToPng(const std::string& path, const Image& image);

}

#endif  // IMAGE_SAVE_TO_PNG_HH_