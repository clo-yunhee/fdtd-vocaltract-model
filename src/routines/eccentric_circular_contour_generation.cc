#include "eccentric_circular_contour_generation.hh"

#include "boundary_interpolation.hh"
#include "tube_segment_connector.hh"

Array2X<uint32_t> vt::eccentricCircularContourGeneration(
    Tensor4<double>& PV_N, uint32_t numSections,
    uint32_t                          totalTubeLengthInCells,
    const Ref<const ArrayX<uint32_t>> tubeSectionDiameterCells,
    const Ref<const ArrayXd> tubeCummSectionLength, bool boundaryInterpolation,
    StartInfo tubeStart, bool pmlSwitch, uint32_t pmlLayer, double ds) {
    // Set a counter to traverse through tubeCummSectionLength
    uint32_t     tubeSegmentCounter = 1;
    const double cellHalfLen = ds / 2;

    // Define currTubeSectionDiameterCells
    auto currTubeSectionDiameterCells_SegmentCounter =
        zeros<Array2X<uint32_t>>(2, totalTubeLengthInCells);
    auto tubeRadiusArray = zeros<ArrayX<uint32_t>>(totalTubeLengthInCells);

    // Set starting coordinate of tube mid-sagittal axis
    // For eccentric geometry we don't need startY
    const uint32_t startX = tubeStart.startX;
    const uint32_t startZ = tubeStart.startZ;

    // Each tube segment consists of number of yz-planes. We first store the
    // tube radius for each of those planes.
    uint32_t currGetRadius;
    uint32_t prevGetRadius = 0;

    for (uint32_t tubeLenCellsCount = 1;
         tubeLenCellsCount <= totalTubeLengthInCells; ++tubeLenCellsCount) {
        // Check the current tube length
        const double currTubeLength = tubeLenCellsCount * ds;

        // Verify if the currTubeLength is more than the tubeCummSectionLength
        // for the current sectionCounter
        // if small or equal then set the tube wall as expected-normal case
        if (currTubeLength <= tubeCummSectionLength(tubeSegmentCounter)) {
            // Get the tube radius
            // We are subtracting 1 as we'll assume that there is a middle
            // row of cells which will act like a mirror/centerline.
            currGetRadius =
                (tubeSectionDiameterCells(tubeSegmentCounter) - 1) / 2;

            // Store the current cross-section tube segment counter
            currTubeSectionDiameterCells_SegmentCounter(2, tubeLenCellsCount) =
                tubeSegmentCounter;
        } else {
            // If the currTubeLength is greater than the actual
            // tubeCummSectionLength Find the difference between currTubeLength
            // and tubeCummSectionLength
            const double diffLength =
                currTubeLength - tubeCummSectionLength(tubeSegmentCounter);

            if (diffLength > cellHalfLen && tubeSegmentCounter != numSections) {
                // Increase the tubeSegmentCounter
                tubeSegmentCounter++;

                // Get the radius for that cross-section
                // [-1: each segment consists of an odd number of cells]
                currGetRadius =
                    (tubeSectionDiameterCells(tubeSegmentCounter) - 1) / 2;

                // Store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(
                    2, tubeLenCellsCount) = tubeSegmentCounter;
            } else {
                // Get the radius for that cross-section
                // [-1: each segment consists of an odd number of cells]
                currGetRadius =
                    (tubeSectionDiameterCells(tubeSegmentCounter) - 1) / 2;

                // Store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(
                    2, tubeLenCellsCount) = tubeSegmentCounter;

                // And then increase the tubeSegmentCounter
                tubeSegmentCounter++;
            }
        }

        tubeRadiusArray(tubeLenCellsCount) = currGetRadius;
    }

    if (boundaryInterpolation) {
        tubeRadiusArray = vt::boundaryInterpolation(tubeRadiusArray, ds);
    }

    // Store the current cross-section tube segment diameter
    currTubeSectionDiameterCells_SegmentCounter.row(1) =
        tubeRadiusArray * 2 + 1;

    // Now as we have tubeRadius for each yz-plane, draw the vocal tract
    // contour starting from the tubeStart.startX

    for (uint32_t tubeLenCellsCount = 1;
         tubeLenCellsCount <= totalTubeLengthInCells; ++tubeLenCellsCount) {
        const uint32_t currGetRadius = tubeRadiusArray(tubeLenCellsCount);

        // For eccentric configuration, the y coordinate of the tube central
        // axis will change but the z coordinate will remain the same for
        // all the tube segments
        // And the X axis will change as we move away from the
        // glottal-end. We need to draw a circle in the yz-plane

        // [-1 to start constructing tube exactly from the tubeStart.startX
        // position]

        const uint32_t tubeX = startX + (tubeLenCellsCount - 1);
        // [for each iteration tubeX remains fixed]

        uint32_t tubeUpperY =
            1 + (pmlLayer * pmlSwitch) + 1;  // deadcell+pmllayers+tubewall
        uint32_t tubeLowerY = tubeUpperY + (currGetRadius * 2 + 1) +
                              1;  // tubeUpperY+tubeDiameter+nextcell
        uint32_t tubeMidY = tubeUpperY + currGetRadius + 1;

        // tubeRightZ and tubeLeftZ mean the right and left side of
        // the center in yz-plane.
        uint32_t tubeRightZ = startZ;
        uint32_t tubeLeftZ = startZ;

        PV_N(tubeUpperY, tubeX, tubeRightZ, 5) = vt::cell_wall;
        PV_N(tubeLowerY, tubeX, tubeLeftZ, 5) = vt::cell_wall;

        // Store the yCoordinate to close if there are any gaps
        uint32_t prevTubeUpperY = tubeUpperY;
        uint32_t prevTubeLowerY = tubeLowerY;

        for (uint32_t radiusCounter = 1; radiusCounter <= currGetRadius;
             ++radiusCounter) {
            const double rActual = (currGetRadius + 0.5) * ds;
            const double zPrime = radiusCounter * ds;

            const double currHeightFromCentralAxis =
                sqrt(rActual * rActual - zPrime * zPrime);
            const double currHeightFromTopOfCentralAxisRowCell =
                currHeightFromCentralAxis - (ds / 2);

            const uint32_t currHeightInCells =
                std::round(currHeightFromTopOfCentralAxisRowCell / ds);

            // Define tube cross-section coordinates in the yz plane
            tubeRightZ = startZ + radiusCounter;
            tubeLeftZ = startZ - radiusCounter;

            tubeUpperY = tubeMidY - currHeightInCells - 1;
            tubeLowerY = tubeMidY + currHeightInCells + 1;

            PV_N(tubeUpperY, tubeX, tubeRightZ, 5) = vt::cell_wall;
            PV_N(tubeLowerY, tubeX, tubeRightZ, 5) = vt::cell_wall;

            PV_N(tubeUpperY, tubeX, tubeLeftZ, 5) = vt::cell_wall;
            PV_N(tubeLowerY, tubeX, tubeLeftZ, 5) = vt::cell_wall;

            // Connect grid cells in yz-plane if there are any gaps
            if (tubeUpperY - prevTubeUpperY > 1) {
                for (uint32_t y = prevTubeUpperY + 1; y <= tubeUpperY; ++y) {
                    PV_N(y, tubeX, tubeRightZ, 5) = vt::cell_wall;
                    PV_N(y, tubeX, tubeLeftZ, 5) = vt::cell_wall;
                }
                for (uint32_t y = tubeLowerY; y <= prevTubeLowerY - 1; ++y) {
                    PV_N(y, tubeX, tubeRightZ, 5) = vt::cell_wall;
                    PV_N(y, tubeX, tubeLeftZ, 5) = vt::cell_wall;
                }
            }

            prevTubeUpperY = tubeUpperY;
            prevTubeLowerY = tubeLowerY;
        }

        // In the end close the circle by using walls at the complete right
        // and left side of the yz-plane [as per the tube radius]

        tubeUpperY = tubeMidY;
        tubeLowerY = tubeMidY;
        tubeRightZ = startZ + currGetRadius + 1;
        tubeLeftZ = startZ - currGetRadius - 1;

        PV_N(tubeUpperY, tubeX, tubeRightZ, 5) = vt::cell_wall;
        PV_N(tubeLowerY, tubeX, tubeLeftZ, 5) = vt::cell_wall;

        if (tubeUpperY - prevTubeUpperY > 1) {
            for (uint32_t y = prevTubeUpperY + 1; y <= tubeUpperY; ++y) {
                PV_N(y, tubeX, tubeRightZ, 5) = vt::cell_wall;
                PV_N(y, tubeX, tubeLeftZ, 5) = vt::cell_wall;
            }
            for (uint32_t y = tubeLowerY; y <= prevTubeLowerY - 1; ++y) {
                PV_N(y, tubeX, tubeRightZ, 5) = vt::cell_wall;
                PV_N(y, tubeX, tubeLeftZ, 5) = vt::cell_wall;
            }
        }

        // Connect grid cells between two adjacent yz-planes if there are any
        // gaps
        if (tubeX != startX && prevGetRadius != currGetRadius) {
            vt::tubeSegmentConnector(PV_N, prevGetRadius, currGetRadius, tubeX);
        }

        // Store the grid radius for the current yz-plane
        prevGetRadius = currGetRadius;
    }

    return currTubeSectionDiameterCells_SegmentCounter;
}