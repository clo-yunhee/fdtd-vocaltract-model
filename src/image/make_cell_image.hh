#ifndef IMAGE_MAKE_CELL_IMAGE_HH_
#define IMAGE_MAKE_CELL_IMAGE_HH_

#include "../types.hh"
#include "image.hh"

namespace image {

Image makeCellImage(const Tensor2<vt::GridCellType>& frame);

}  // namespace image

#endif  // IMAGE_MAKE_CELL_IMAGE_HH_