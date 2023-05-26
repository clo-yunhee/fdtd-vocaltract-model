#include "eccentric_elliptical_contour_generation.hh"

#include "boundary_interpolation.hh"
#include "tube_segment_connector.hh"

array vt::eccentricEllipticalContourGeneration(
    array& PV_N, uint32_t numSections, uint32_t totalTubeLengthInCells,
    const array& tubeCummSectionLength, const array& ellipseAxisLenInfo,
    bool boundaryInterpolation, StartInfo tubeStart, bool pmlSwitch,
    uint32_t pmlLayer, double ds) {
    // Set a counter to traverse through tubeCummSectionLength
    uint32_t     tubeSegmentCounter = 0;
    const double cellHalfLen = ds / 2;

    // Define currTubeSectionDiameterCells
    array currTubeSectionDiameterCells_SegmentCounter =
        constant(0, 2, totalTubeLengthInCells);
    array tubeSemiMinorAxisRadiusInCells =
        constant(0, 1, totalTubeLengthInCells, u32);
    array tubeSemiMajorAxisRadiusInCells =
        constant(0, 1, totalTubeLengthInCells, u32);

    // Set starting coordinate of tube mid-sagittal axis
    // For eccentric geometry we don't need startY
    const uint32_t startX = tubeStart.startX;
    const uint32_t startZ = tubeStart.startZ;

    // Each tube segment consists of number of yz-planes. We first store the
    // tube radius for each of those planes.
    uint32_t currSemiMinorAxisRadius;
    uint32_t currSemiMajorAxisRadius;
    uint32_t prevSemiMajorAxisRadius = 0;

    for (uint32_t tubeLenCellsCount = 0;
         tubeLenCellsCount < totalTubeLengthInCells; ++tubeLenCellsCount) {
        // Check the current tube length
        const double currTubeLength = (tubeLenCellsCount + 1) * ds;

        // Verify if the currTubeLength is more than the tubeCummSectionLength
        // for the current sectionCounter
        // if small or equal then set the tube wall as expected-normal case
        if (allTrue<bool>(currTubeLength <=
                          tubeCummSectionLength(tubeSegmentCounter))) {
            // Get the tube radius
            // We are subtracting 1 as we'll assume that there is a middle
            // row of cells which will act like a mirror/centerline.
            currSemiMajorAxisRadius =
                ((ellipseAxisLenInfo(0, tubeSegmentCounter) - 1) / 2)
                    .scalar<uint32_t>();
            currSemiMinorAxisRadius =
                ((ellipseAxisLenInfo(1, tubeSegmentCounter) - 1) / 2)
                    .scalar<uint32_t>();

            // Store the current cross-section tube segment counter
            currTubeSectionDiameterCells_SegmentCounter(1, tubeLenCellsCount) =
                tubeSegmentCounter;
        } else {
            // If the currTubeLength is greater than the actual
            // tubeCummSectionLength Find the difference between currTubeLength
            // and tubeCummSectionLength
            const double diffLength =
                currTubeLength -
                tubeCummSectionLength(tubeSegmentCounter).scalar<float>();

            if (diffLength > cellHalfLen && tubeSegmentCounter != numSections) {
                // Increase the tubeSegmentCounter
                tubeSegmentCounter = tubeSegmentCounter + 1;

                // Get the radius for that cross-section
                // [-1: each segment consists of an odd number of cells]
                currSemiMajorAxisRadius =
                    ((ellipseAxisLenInfo(0, tubeSegmentCounter) - 1) / 2)
                        .scalar<uint32_t>();
                currSemiMinorAxisRadius =
                    ((ellipseAxisLenInfo(1, tubeSegmentCounter) - 1) / 2)
                        .scalar<uint32_t>();

                // Store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(
                    1, tubeLenCellsCount) = tubeSegmentCounter;
            } else {
                // Get the radius for that cross-section
                // [-1: each segment consists of an odd number of cells]
                currSemiMajorAxisRadius =
                    ((ellipseAxisLenInfo(0, tubeSegmentCounter) - 1) / 2)
                        .scalar<uint32_t>();
                currSemiMinorAxisRadius =
                    ((ellipseAxisLenInfo(1, tubeSegmentCounter) - 1) / 2)
                        .scalar<uint32_t>();

                // Store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(
                    1, tubeLenCellsCount) = tubeSegmentCounter;

                // And then increase the tubeSegmentCounter
                tubeSegmentCounter = tubeSegmentCounter + 1;
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
    currTubeSectionDiameterCells_SegmentCounter(0, span) =
        tubeSemiMinorAxisRadiusInCells * 2 + 1;

    // Now as we have tubeRadius for each yz-plane, draw the vocal tract
    // contour starting from the tubeStart.startX

    for (uint32_t tubeLenCellsCount = 0;
         tubeLenCellsCount < totalTubeLengthInCells; ++tubeLenCellsCount) {
        const uint32_t currSemiMajorAxisRadius =
            tubeSemiMajorAxisRadiusInCells(tubeLenCellsCount)
                .scalar<uint32_t>();
        const uint32_t currSemiMinorAxisRadius =
            tubeSemiMinorAxisRadiusInCells(tubeLenCellsCount)
                .scalar<uint32_t>();

        // The Y and Z coordinate of the tube central axis is not gonna
        // change. Only the X axis will change as we move away from the
        // glottal-end. We need to draw an ellipse in the yz-plane

        const uint32_t tubeX = startX + tubeLenCellsCount;
        // [for each iteration tubeX remains fixed]

        uint32_t tubeUpperY =
            1 + (pmlLayer * pmlSwitch) + 1;  // deadcell+pmllayers+tubewall
        uint32_t tubeLowerY = tubeUpperY + (currSemiMinorAxisRadius * 2 + 1) +
                              1;  // tubeUpperY+tubeDiameter+nextcell
        uint32_t tubeMidY = tubeUpperY + currSemiMinorAxisRadius + 1;

        // tubeRightZ and tubeLeftZ mean the right and left side of
        // the center in yz-plane.
        uint32_t tubeRightZ = startZ;
        uint32_t tubeLeftZ = startZ;

        PV_N(tubeUpperY, tubeX, tubeRightZ, 4) = vt::cell_wall;
        PV_N(tubeLowerY, tubeX, tubeLeftZ, 4) = vt::cell_wall;

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

            // Define tubeUpperY and tubeLowerY
            tubeUpperY = tubeMidY - currHeightInCells - 1;
            tubeLowerY = tubeMidY + currHeightInCells + 1;

            PV_N(tubeUpperY, tubeX, tubeRightZ, 4) = vt::cell_wall;
            PV_N(tubeLowerY, tubeX, tubeRightZ, 4) = vt::cell_wall;

            PV_N(tubeUpperY, tubeX, tubeLeftZ, 4) = vt::cell_wall;
            PV_N(tubeLowerY, tubeX, tubeLeftZ, 4) = vt::cell_wall;

            // Connect grid cells in yz-plane if there are any gaps
            if (tubeUpperY - prevTubeUpperY > 1) {
                PV_N(seq(prevTubeUpperY + 1, tubeUpperY), tubeX, tubeRightZ,
                     4) = vt::cell_wall;
                PV_N(seq(prevTubeUpperY + 1, tubeUpperY), tubeX, tubeLeftZ, 4) =
                    vt::cell_wall;
                PV_N(seq(tubeLowerY, prevTubeLowerY - 1), tubeX, tubeRightZ,
                     4) = vt::cell_wall;
                PV_N(seq(tubeLowerY, prevTubeLowerY - 1), tubeX, tubeLeftZ, 4) =
                    vt::cell_wall;
            }

            prevTubeUpperY = tubeUpperY;
            prevTubeLowerY = tubeLowerY;
        }

        // In the end close the circle by using walls at the complete right
        // and left side of the yz-plane [as per the tube radius]

        tubeUpperY = tubeMidY;
        tubeLowerY = tubeMidY;
        tubeRightZ = startZ + currSemiMajorAxisRadius + 1;
        tubeLeftZ = startZ - currSemiMajorAxisRadius - 1;

        PV_N(tubeUpperY, tubeX, tubeRightZ, 4) = vt::cell_wall;
        PV_N(tubeLowerY, tubeX, tubeLeftZ, 4) = vt::cell_wall;

        if (tubeUpperY - prevTubeUpperY > 1) {
            PV_N(seq(prevTubeUpperY + 1, tubeUpperY), tubeX, tubeRightZ, 4) =
                vt::cell_wall;
            PV_N(seq(prevTubeUpperY + 1, tubeUpperY), tubeX, tubeLeftZ, 4) =
                vt::cell_wall;
            PV_N(seq(tubeLowerY, prevTubeLowerY - 1), tubeX, tubeRightZ, 4) =
                vt::cell_wall;
            PV_N(seq(tubeLowerY, prevTubeLowerY - 1), tubeX, tubeLeftZ, 4) =
                vt::cell_wall;
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