#include "make_wave_image.hh"

#include <colormap/palettes.hpp>

using namespace vt;

image::Image image::makeWaveImage(const Tensor2<double>& frame,
                                  const double min, const double max) {
    const auto frameSize = dimensions(frame);

    const uint32_t frameWidth = frameSize[2];
    const uint32_t frameHeight = frameSize[1];

    const uint32_t cellSize = 8;

    const uint32_t width = frameWidth * cellSize;
    const uint32_t height = frameHeight * cellSize;

    image::Image image{};
    image.width = width;
    image.height = height;
    image.pixels.resize(height * width);

    const auto map = colormap::palettes.at("inferno");

    for (uint32_t fy = 0; fy < frameHeight; ++fy) {
        for (uint32_t fx = 0; fx < frameWidth; ++fx) {
            double val = frame(fy + 1, fx + 1);
            if (val < min) val = min;
            if (val > max) val = max;

            const double brightness = ((val - min) / (max - min));
            const auto&  color = map(brightness);

            const uint32_t x0 = fx * cellSize;
            const uint32_t y0 = fy * cellSize;

            for (uint32_t y = 0; y < cellSize; ++y) {
                for (uint32_t x = 0; x < cellSize; ++x) {
                    image.pixels[(y0 + y) * width + (x0 + x)].r =
                        color.getRed().getValue();
                    image.pixels[(y0 + y) * width + (x0 + x)].g =
                        color.getGreen().getValue();
                    image.pixels[(y0 + y) * width + (x0 + x)].b =
                        color.getBlue().getValue();
                    image.pixels[(y0 + y) * width + (x0 + x)].a = 0xFF;
                }
            }
        }
    }

    return image;
}
