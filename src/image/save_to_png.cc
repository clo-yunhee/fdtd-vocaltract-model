#include "save_to_png.hh"

#include <png.h>

#include <cstdio>
#include <memory>

void image::saveToPng(const std::string& path, const Image& image) {
    const uint32_t width = image.width;
    const uint32_t height = image.height;

    FILE* fp = fopen(path.c_str(), "wb");
    if (!fp) {
        perror("png fopen failed:");
        return;
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr,
                                              nullptr, nullptr);
    if (!png) {
        fprintf(stderr, "png_create_write_struct failed\n");
        return;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "png_create_info_struct failed\n");
        return;
    }

    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "png setjmp failed\n");
        return;
    }

    png_init_io(png, fp);

    // Output is 8bit depth, RGBA format.
    png_set_IHDR(png, info, image.width, image.height, 8, PNG_COLOR_TYPE_RGBA,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    std::vector<std::unique_ptr<png_byte[]>> row_pointers_u(height);
    std::unique_ptr<png_bytep[]>             row_pointers =
        std::make_unique<png_bytep[]>(height);

    for (uint32_t y = 0; y < height; ++y) {
        row_pointers_u[y] = std::make_unique<png_byte[]>(width * 4);
        row_pointers[y] = row_pointers_u[y].get();

        png_bytep row = row_pointers[y];
        for (uint32_t x = 0; x < width; ++x) {
            png_bytep px = &(row[x * 4]);
            px[0] = image.pixels[y * width + x].r;
            px[1] = image.pixels[y * width + x].g;
            px[2] = image.pixels[y * width + x].b;
            px[3] = image.pixels[y * width + x].a;
        }
    }

    png_write_image(png, row_pointers.get());
    png_write_end(png, nullptr);

    fclose(fp);

    png_destroy_write_struct(&png, &info);
}