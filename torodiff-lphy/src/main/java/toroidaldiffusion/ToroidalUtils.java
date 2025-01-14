package toroidaldiffusion;

public class ToroidalUtils {
    public static double wrapToMaxAngle(double rawAngle, double MAX_ANGLE_VALUE) {
        if (rawAngle > MAX_ANGLE_VALUE) {
            int K = (int)Math.floor(rawAngle / MAX_ANGLE_VALUE);
            double fractionRemainder = rawAngle / MAX_ANGLE_VALUE - K;
            return fractionRemainder * MAX_ANGLE_VALUE;
        }

        if (rawAngle < 0.0) {
            int K = (int)Math.floor(-rawAngle / MAX_ANGLE_VALUE);
            double fractionRemainder = (-rawAngle / MAX_ANGLE_VALUE) - K;
            return MAX_ANGLE_VALUE - (fractionRemainder * MAX_ANGLE_VALUE);
        }

        return rawAngle;
    }

    /**
     * This replicate the %% operator in R,
     * which gives the same result as {@link #wrapToMaxAngle(double, double)},
     * when lower = 0.
     * @param x   dividend, e.g. rawAngle
     * @param y   divisor, e.g. MAX_ANGLE_VALUE
     * @return    The result of %% operator in R.
     *            It is useful for operations where wrap-around behavior is needed.
     */
    public static double moduloInR(double x, double y) {
        return ((x % y) + Math.abs(y)) % Math.abs(y);
    }

}
