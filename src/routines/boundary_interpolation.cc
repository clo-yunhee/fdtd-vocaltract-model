#include "boundary_interpolation.hh"

array vt::boundaryInterpolation(const array& tubeRadiusArray, double ds) {
    // Illustration:
    //                     C         C
    //                     *         *
    //             E   *   *         *   *   E
    // 		  	  *	  *   *   or    *   *   *
    // 		  *	  *	  *   *         *   *   *   *
    // 	  *	  *	  *	  *   *         *   *   *   *   *
    // *********************         *********************
    // A           D       B         B       D           A
    //         (Case 1)         or          (case 2)

    // In the Triangle ABC, let's assume, we need to find DE = ?
    // AD/AB = DE/BC => DE = (AD*BC)/AB; [Similar Triangle Equality]

    // Create a new tubeRadiusArray to store the radius of the interpolated
    // boundary.
    array tubeNewRadiusArray = constant(0, 1, tubeRadiusArray.elements());

    // Set a counter to point to the start of the triangle
    uint32_t triangleStartCounter = 0;

    while (triangleStartCounter + 1 < tubeRadiusArray.elements()) {
        // Set a counter and traverse from start of the tube to
        // the end of the tube
        uint32_t tubeCounter = 1;

        // To use the "Similar Triangle Equality" first find
        // the start and end of the triangle
        while (allTrue<bool>(
            tubeRadiusArray(triangleStartCounter) ==
            tubeRadiusArray(triangleStartCounter + tubeCounter))) {
            if (triangleStartCounter + tubeCounter !=
                tubeRadiusArray.elements()) {
                tubeCounter = tubeCounter + 1;
            } else {
                break;
            }
        }

        // Set triangle length and height
        const uint32_t triangleLengthInCells = tubeCounter;
        const int32_t  triangleHeightInCells =
            (tubeRadiusArray(triangleStartCounter + tubeCounter)
                 .scalar<int32_t>() -
             tubeRadiusArray(triangleStartCounter).scalar<int32_t>());

        // Set the tubeNewRadiusArray for the triangleStart and triangleEnd
        tubeNewRadiusArray(triangleStartCounter) =
            tubeRadiusArray(triangleStartCounter);
        tubeNewRadiusArray(triangleStartCounter + tubeCounter) =
            tubeRadiusArray(triangleStartCounter + tubeCounter);

        if (triangleHeightInCells > 0) {
            // Case 1 - Check the illustration above
            const double AB = triangleLengthInCells * ds;
            const double BC = triangleHeightInCells * ds;

            for (uint32_t heightCounter = 1; heightCounter <= tubeCounter - 1;
                 ++heightCounter) {
                const double   AD = heightCounter * ds;
                const double   DE = (AD * BC) / AB;
                const uint32_t numCells = (uint32_t)std::round(DE / ds);

                tubeNewRadiusArray(triangleStartCounter + heightCounter) =
                    tubeRadiusArray(triangleStartCounter + heightCounter) +
                    numCells;
            }
        } else {
            // Case 2 - Check the illustration above
            const double AB = triangleLengthInCells * ds;
            const double BC = std::abs(triangleHeightInCells) * ds;

            for (uint32_t heightCounter = 1; heightCounter <= tubeCounter - 1;
                 ++heightCounter) {
                const double   AD = (tubeCounter - heightCounter) * ds;
                const double   DE = (AD * BC) / AB;
                const uint32_t numCells = (uint32_t)std::round(DE / ds);

                tubeNewRadiusArray(triangleStartCounter + heightCounter) =
                    tubeRadiusArray(triangleStartCounter + heightCounter) +
                    numCells;
            }
        }

        // Set the triangleStartCounter
        triangleStartCounter += tubeCounter;
    }

    return tubeNewRadiusArray;
}