#include "find_cell_types.hh"

array vt::findCellTypes(const array& PV_N, const uint32_t tubeX) {
    // For a given yz-plane, this function determines whether a particular
    // grid cell is inside the tube contour, outside of the tube contour or
    // on the tube contour [tubeContour].

    // Define grid cells as inside the VT contour, outside the VT contour or
    // on the VT contour.

    const dim4 gridSize = PV_N.dims();

    array gridPlaneProp = constant(0, gridSize[0], gridSize[2]);

    for (uint32_t zCount = 0; zCount < gridSize[2]; ++zCount) {
        const array findPlaneWalls =
            where(PV_N(span, tubeX, zCount, 4) == (uint32_t)vt::cell_wall);

        gridPlaneProp(findPlaneWalls, zCount) =
            (float)GridCellTypeInplane::onVTContour;

        if (findPlaneWalls.isempty()) {
            gridPlaneProp(span, zCount) =
                (float)GridCellTypeInplane::outVTContour;
        } else {
            const auto minWall = min<uint32_t>(findPlaneWalls);
            const auto maxWall = max<uint32_t>(findPlaneWalls);

            gridPlaneProp(seq(0, minWall - 1), zCount) =
                (float)GridCellTypeInplane::outVTContour;
            if (maxWall + 1 <= gridSize[0] - 1)
                gridPlaneProp(seq(maxWall + 1, gridSize[0] - 1), zCount) =
                    (float)GridCellTypeInplane::outVTContour;

            const array inVTContourCells = gridPlaneProp(span, zCount) == 0;
            gridPlaneProp(inVTContourCells, zCount) =
                (float)GridCellTypeInplane::inVTContour;
        }
    }

    return gridPlaneProp;
}