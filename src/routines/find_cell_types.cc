#include "find_cell_types.hh"

ArrayXX<vt::GridCellTypeInplane> vt::findCellTypes(const Tensor4<double>& PV_N,
                                                   const uint32_t tubeX) {
    // For a given yz-plane, this function determines whether a particular
    // grid cell is inside the tube contour, outside of the tube contour or
    // on the tube contour [tubeContour].

    // Define grid cells as inside the VT contour, outside the VT contour or
    // on the VT contour.

    const auto gridSize = dimensions(PV_N);
    auto       gridPlaneProp =
        zeros<ArrayXX<GridCellTypeInplane>>(gridSize[1], gridSize[3]);

    gridPlaneProp.setConstant(GridCellTypeInplane::null);

    for (uint32_t zCount = 1; zCount <= gridSize[3]; ++zCount) {
        // Find plane walls
        bool     planeWallsFound(false);
        uint32_t minWall(UINT32_MAX);
        uint32_t maxWall(0);

        for (uint32_t y = 1; y <= gridSize[1]; ++y) {
            if (PV_N(y, tubeX, zCount, 5) == vt::cell_wall) {
                gridPlaneProp(y, zCount) = GridCellTypeInplane::onVTContour;

                if (!planeWallsFound) planeWallsFound = true;
                if (y < minWall) minWall = y;
                if (y > maxWall) maxWall = y;
            }
        }

        if (!planeWallsFound) {
            gridPlaneProp(seq(1, gridSize[1]), zCount) =
                GridCellTypeInplane::outVTContour;
        } else {
            gridPlaneProp(seq(1, minWall - 1), zCount) =
                GridCellTypeInplane::outVTContour;
            gridPlaneProp(seq(maxWall + 1, gridSize[1]), zCount) =
                GridCellTypeInplane::outVTContour;

            for (uint32_t y = 1; y <= gridSize[1]; ++y) {
                if (gridPlaneProp(y, zCount) == GridCellTypeInplane::null) {
                    gridPlaneProp(y, zCount) = GridCellTypeInplane::inVTContour;
                }
            }
        }
    }

    return gridPlaneProp;
}