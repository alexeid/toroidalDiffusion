package toroidaldiffusion.operator;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class WrappedRandomWalkOperatorTest {

    @Test
    void proposal() {
    }

    private static final double TOLERANCE = 1e-10;

    @Test
    void wrapAngles() {
        double[] withinBounds = {0.0, Math.PI / 4, -Math.PI / 4, Math.PI, -Math.PI, Math.PI / 2, -Math.PI / 2, 1.0, -1.0, 2.0};
        double[] expectedWithinBounds = {0.0, Math.PI / 4, -Math.PI / 4, Math.PI, -Math.PI, Math.PI / 2, -Math.PI / 2, 1.0, -1.0, 2.0};

        // Angles outside the bounds of -π to π
        double[] outsideBounds = {Math.PI + 0.1, -Math.PI - 0.1, 2 * Math.PI, -2 * Math.PI, 3 * Math.PI - 0.1, -3 * Math.PI + 0.1, 4 * Math.PI, -4 * Math.PI, 5 * Math.PI, -5 * Math.PI};
        double[] expectedOutsideBounds = {-Math.PI + 0.1, Math.PI - 0.1, 0.0, 0.0, Math.PI - 0.1, -Math.PI + 0.1, 0.0, 0.0, Math.PI, -Math.PI};

        // Test angles within bounds
        for (int i = 0; i < withinBounds.length; i++) {
            double observed = WrappedRandomWalkOperator.wrapAngles(withinBounds[i]);
            assertEquals(expectedWithinBounds[i], observed, TOLERANCE, "Failed for angle: " + withinBounds[i]);
        }

        // Test angles outside bounds
        for (int i = 0; i < outsideBounds.length; i++) {
            double observed = WrappedRandomWalkOperator.wrapAngles(outsideBounds[i]);
            assertEquals(expectedOutsideBounds[i], observed, TOLERANCE, "Failed for angle: " + outsideBounds[i]);
        }
    }

    }
