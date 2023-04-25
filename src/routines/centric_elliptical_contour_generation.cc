#include "centric_elliptical_contour_generation.hh"

#include "boundary_interpolation.hh"
#include "tube_segment_connector.hh"

Array2X<uint32_t> vt::centricEllipticalContourGeneration(
    Tensor4<double>& PV_N, uint32_t numSections,
    uint32_t                           totalTubeLengthInCells,
    const Ref<const ArrayXd>           tubeCummSectionLength,
    const Ref<const Array2X<uint32_t>> ellipseAxisLenInfo,
    bool boundaryInterpolation, StartInfo tubeStart, double ds) {
    // Set a counter to traverse through tubeCummSectionLength
    uint32_t     tubeSegmentCounter = 1;
    const double cellHalfLen = ds / 2;

    // Define currTubeSectionDiameterCells
    // For ellptical cross-section, we will use this array to store the minor
    // axis tube diamaters.
    auto currTubeSectionDiameterCells_SegmentCounter =
        zeros<Array2X<uint32_t>>(2, totalTubeLengthInCells);
    auto tubeSemiMinorAxisRadiusInCells =
        zeros<ArrayX<uint32_t>>(totalTubeLengthInCells);
    auto tubeSemiMajorAxisRadiusInCells =
        zeros<ArrayX<uint32_t>>(totalTubeLengthInCells);

    // Set starting coordinate of tube mid-sagittal axis
    const uint32_t startX = tubeStart.startX;
    const uint32_t startY = tubeStart.startY;
    const uint32_t startZ = tubeStart.startZ;

    // Each tube segment consists of number of yz-planes. We first store the
    // tube semi-minor axis radius/legnth for each of those planes.
    uint32_t currSemiMinorAxisRadius;
    uint32_t currSemiMajorAxisRadius;
    uint32_t prevSemiMajorAxisRadius = 0;

    for (uint32_t tubeLenCellsCount = 1;
         tubeLenCellsCount <= totalTubeLengthInCells; ++tubeLenCellsCount) {
        // Check the current tube length
        const double currTubeLength = tubeLenCellsCount * ds;

        // Verify if the currTubeLength is more than the tubeCummSectionLength
        // for the current sectionCounter
        // if small or equal then set the tube wall as expected-normal Case
        if (currTubeLength <= tubeCummSectionLength(tubeSegmentCounter)) {
            // Get the tube radius
            // We are subtracting 1 as we'll assume that there is a middle
            // row of cells which will act like a mirror/centerline.
            currSemiMajorAxisRadius =
                (ellipseAxisLenInfo(1, tubeSegmentCounter) - 1) / 2;
            currSemiMinorAxisRadius =
                (ellipseAxisLenInfo(2, tubeSegmentCounter) - 1) / 2;

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
                currSemiMajorAxisRadius =
                    (ellipseAxisLenInfo(1, tubeSegmentCounter) - 1) / 2;
                currSemiMinorAxisRadius =
                    (ellipseAxisLenInfo(2, tubeSegmentCounter) - 1) / 2;

                // Store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(
                    2, tubeLenCellsCount) = tubeSegmentCounter;
            } else {
                // Get the radius for that cross-section
                // [-1: each segment consists of an odd number of cells]
                currSemiMajorAxisRadius =
                    (ellipseAxisLenInfo(1, tubeSegmentCounter) - 1) / 2;
                currSemiMinorAxisRadius =
                    (ellipseAxisLenInfo(2, tubeSegmentCounter) - 1) / 2;

                // Store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(
                    2, tubeLenCellsCount) = tubeSegmentCounter;

                // And then increase the tubeSegmentCounter
                tubeSegmentCounter++;
            }
        }

        tubeSemiMinorAxisRadiusInCells(tubeLenCellsCount) =
            currSemiMinorAxisRadius;
        tubeSemiMajorAxisRadiusInCells(tubeLenCellsCount) =
            currSemiMajorAxisRadius;
    }

    if (boundaryInterpolation) {
        // TODO: Deb - Think of what's the correct method to interpolate
        // between two tube segmnts.
    }

    // Store the current cross-section tube segment diameter
    currTubeSectionDiameterCells_SegmentCounter.row(1) =
        tubeSemiMinorAxisRadiusInCells * 2 + 1;

    // Now as we have tubeRadius for each yz-plane, draw the vocal tract
    // contour starting from the tubeStart.startX

    for (uint32_t tubeLenCellsCount = 1;
         tubeLenCellsCount <= totalTubeLengthInCells; ++tubeLenCellsCount) {
        const uint32_t currSemiMajorAxisRadius =
            tubeSemiMajorAxisRadiusInCells(tubeLenCellsCount);
        const uint32_t currSemiMinorAxisRadius =
            tubeSemiMinorAxisRadiusInCells(tubeLenCellsCount);

        // The Y and Z coordinate of the tube central axis is not gonna
        // change. Only the X axis will change as we move away from the
        // glottal-end. We need to draw an ellipse in the yz-plane

        // [-1 to start constructing tube exactly from the tubeStart.startX
        // position]

        const uint32_t tubeX = startX + (tubeLenCellsCount - 1);
        // [for each iteration tubeX remains fixed]

        uint32_t tubeUpperY = startY - currSemiMinorAxisRadius - 1;
        uint32_t tubeLowerY = startY + currSemiMinorAxisRadius + 1;

        // tubeRightZ and tubeLeftZ mean the right and left side of
        // the center in yz-plane.
        uint32_t tubeRightZ = startZ;
        uint32_t tubeLeftZ = startZ;

        PV_N(tubeUpperY, tubeX, tubeRightZ, 5) = vt::cell_wall;
        PV_N(tubeLowerY, tubeX, tubeLeftZ, 5) = vt::cell_wall;

        // Store the yCoordinate to close if there are any gaps
        uint32_t prevTubeUpperY = tubeUpperY;
        uint32_t prevTubeLowerY = tubeLowerY;

        for (uint32_t semiMajorAxisRadiusCounter = 1;
             semiMajorAxisRadiusCounter <= currSemiMajorAxisRadius;
             ++semiMajorAxisRadiusCounter) {
            /*
            Ellipse equation: z^2/a^2 + y^2/b^2 = 1 [Note: I have mentioned
            z^2 becasuse the semi-major axis is along the z-axis]
            a = (currSemiMajorAxisRadius+0.5)*ds [Length of semiMaorAxisRadius]
            b = (currSemiMinorAxisRadius+0.5)*ds [Length of semiMinorAxisRadius]
            z = semiMajorAxisRadiusCounter*ds
            Determine y (ellipse_height) = ? using Ellipse equation
            [i.e., height of the ellipse along y-axis]
            */

            const double a = (currSemiMajorAxisRadius + 0.5) * ds;
            const double b = (currSemiMinorAxisRadius + 0.5) * ds;
            const double z = semiMajorAxisRadiusCounter * ds;

            const double currElipseHeightFromCentralAxis =
                b * sqrt(1 - (z * z) / (a * a));
            const double currHeightFromTopOfCentralAxisRowCell =
                currElipseHeightFromCentralAxis - (ds / 2);

            const uint32_t currHeightInCells =
                std::round(currHeightFromTopOfCentralAxisRowCell / ds);

            // Define tube cross-section coordinates in the yz plane
            tubeRightZ = startZ + semiMajorAxisRadiusCounter;
            tubeLeftZ = startZ - semiMajorAxisRadiusCounter;

            tubeUpperY = startY - currHeightInCells - 1;
            tubeLowerY = startY + currHeightInCells + 1;

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

        tubeUpperY = startY;
        tubeLowerY = startY;
        tubeRightZ = startZ + currSemiMajorAxisRadius + 1;
        tubeLeftZ = startZ - currSemiMajorAxisRadius - 1;

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
        if (tubeX != startX &&
            prevSemiMajorAxisRadius != currSemiMajorAxisRadius) {
            vt::tubeSegmentConnector(PV_N, prevSemiMajorAxisRadius,
                                     currSemiMajorAxisRadius, tubeX);
        }

        // Store the grid radius for the current yz-plane
        prevSemiMajorAxisRadius = currSemiMajorAxisRadius;
    }

    return currTubeSectionDiameterCells_SegmentCounter;
}