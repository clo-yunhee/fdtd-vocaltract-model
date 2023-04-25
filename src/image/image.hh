#ifndef IMAGE_IMAGE_HH_
#define IMAGE_IMAGE_HH_

#include <cstdint>
#include <vector>

namespace image {

union Pixel {
    struct {
        uint8_t r;
        uint8_t g;
        uint8_t b;
        uint8_t a;
    };
    uint32_t rgba;
};

struct Image {
    uint32_t           width;
    uint32_t           height;
    std::vector<Pixel> pixels;
};

}  // namespace image

#endif  // IMAGE_IMAGE_HH_