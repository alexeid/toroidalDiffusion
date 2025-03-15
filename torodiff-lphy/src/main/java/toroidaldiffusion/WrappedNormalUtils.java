package toroidaldiffusion;

import lphy.core.model.Value;
import lphy.core.model.ValueUtils;

import static java.lang.Math.sqrt;

public class WrappedNormalUtils {

    /**
     * Reparametrisation: alpha3^2 < alpha1 * alpha2, alpha3 < sqrt(alpha1 * alpha2)
     * @param drift
     * @param driftCorr
     * @return
     */
    public static double[] getAlphaArr(Value<Number[]> drift, Value<Number> driftCorr) {
        double[] twoDrifts = ValueUtils.doubleArrayValue(drift);
        double corr = ValueUtils.doubleValue(driftCorr);

        double alpha1 = twoDrifts[0];
        double alpha2 = twoDrifts[1];
        double alpha3 = sqrt(alpha1*alpha2) * corr;

        return new double[]{alpha1, alpha2, alpha3};
    }

    public static double getDriftCorr(double alpha1, double alpha2, double alpha3) {
        return alpha3 / sqrt(alpha1*alpha2);
    }
}
