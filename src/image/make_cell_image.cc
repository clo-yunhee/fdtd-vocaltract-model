#include "make_cell_image.hh"

using namespace vt;

image::Image image::makeCellImage(const Tensor2<GridCellType>& frame) {
    const auto frameSize = dimensions(frame);

    const uint32_t frameWidth = frameSize[2];
    const uint32_t frameHeight = frameSize[1];

    const uint32_t cellSize = 7;

    const uint32_t width = frameWidth * cellSize;
    const uint32_t height = frameHeight * cellSize;

    image::Image image{};
    image.width = width;
    image.height = height;
    image.pixels.resize(height * width);

    uint32_t rgba;

    for (uint32_t fy = 0; fy < frameHeight; ++fy) {
        for (uint32_t fx = 0; fx < frameWidth; ++fx) {
            const auto cellType = frame(fy + 1, fx + 1);
            switch (cellType) {
                case cell_wall:
                    rgba = 0xFFD700FF;  // yellow-gold
                    break;
                case cell_air:
                    rgba = 0x00000000;  // transparent
                    break;
                case cell_excitation:
                    rgba = 0xFF0000FF;  // red
                    break;
                case cell_dead:
                    rgba = 0x000000FF;  // black
                    break;
                case cell_noPressure:
                    rgba = 0x5A5A5AFF;  // dark gray
                    break;
                default:
                    rgba = 0xFF69B4FF;  // fuschia
                    break;
            }

            const uint32_t x0 = fx * cellSize;
            const uint32_t y0 = fy * cellSize;

            for (uint32_t y = 0; y < cellSize; ++y) {
                for (uint32_t x = 0; x < cellSize; ++x) {
                    image.pixels[(y0 + y) * width + (x0 + x)].rgba = rgba;
                }
            }
        }
    }

    return image;
}
