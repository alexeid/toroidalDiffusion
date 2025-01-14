package toroidaldiffusion;

import lphy.base.distribution.ParametricDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Stream;

import static java.lang.Math.sqrt;
import static toroidaldiffusion.PhyloWrappedBivariateDiffusion.*;

public class WrappedBivariateNormal extends ParametricDistribution<Double[]> {

    Value<Number[]> mu;
    Value<Number[]> sigma;
    Value<Double[]> drift;
    Value<Number> driftCorr;

    NormalDistribution nD1;
    NormalDistribution nD2;

    public WrappedBivariateNormal(@ParameterInfo(name = muParamName, description = "the mean of the stationary distribution.")
                                  Value<Number[]> mu,
                                  @ParameterInfo(name = sigmaParamName, description = "the two variance terms.")
                                  Value<Number[]> sigma,
                                  @ParameterInfo(name = DRIFT_PARAM, description = "the two drift terms : alpha1 and alpha2.")
                                  Value<Double[]> drift,
                                  @ParameterInfo(name = DRIFT_CORR_PARAM, description = "the correlation of two drift terms, " +
                                          "ranged from -1 to 1 excluding -1 and 1, but not alpha3 where alpha3 = sqrt(alpha1*alpha2) * corr.")
                                  Value<Number> driftCorr) {
        super();
        this.mu = mu;
        this.sigma = sigma;
        this.drift = drift;
        this.driftCorr = driftCorr;

        constructDistribution(random);
    }


    @Override
    protected void constructDistribution(RandomGenerator random) {
        // X1 <- rnorm(n, mu1, s1)
        double mu1 = ValueUtils.doubleArrayValue(mu)[0];
        double s1 = ValueUtils.doubleArrayValue(sigma)[0];
        nD1 = new NormalDistribution(random, mu1, s1, NormalDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
    }

    @Override
    public RandomVariable<Double[]> sample() {
        double x1 = nD1.sample();

        double mu1 = ValueUtils.doubleArrayValue(mu)[0];
        double mu2 = ValueUtils.doubleArrayValue(mu)[1];
        // transfer alpha1,2,3 into rho
        double[] alpha = Stream.of(getAlphaarr(drift, driftCorr)).mapToDouble(Double::doubleValue).toArray();
        double rho = getRho(alpha);
        // transfer sigma1,2 into s1,2 for bivariate normal
        double s1 = getS1(ValueUtils.doubleArrayValue(sigma)[0], alpha);
        double s2 = getS2(ValueUtils.doubleArrayValue(sigma)[1], alpha);

        // sample from the second Normal conditioned on the first Normal
        // X2 <- rnorm(n, mu2 + (s2/s1) * rho * (X1 - mu1), sqrt((1 - rho^2)*s2^2))
        double[] mean_sd = getConditionedMeanSd(mu1, s1, mu2, s2, rho, x1);
        nD2 = new NormalDistribution(random, mean_sd[0], mean_sd[1], NormalDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
        double x2 = nD2.sample();

        return new RandomVariable<>("bivariate", new Double[]{x1, x2}, this);
    }

    /**
     * @param alpha  drift term : alpha1, alpha2, alpha3
     * @return  the constant term outside the matrix: 1 / 2(alpha1*alpha2-alpha3*alpha3)
     */
    public static double getConst(double[] alpha) {
        return 0.5 / (alpha[0]*alpha[1]-alpha[2]*alpha[2]);
    }

    public static double getRho(double[] alpha) {
        double constant = getConst(alpha);
        return -alpha[2] * constant / ( sqrt(constant * alpha[1] ) * sqrt(constant * alpha[0] ) );
    }

    public static double getS1(double sigma1, double[] alpha) {
        double constant = getConst(alpha);
        // alpha2 is alpha[1]
        return sqrt(constant * alpha[1]) * sigma1;
    }

    public static double getS2(double sigma2, double[] alpha) {
        double constant = getConst(alpha);
        // alpha1 is alpha[0]
        return sqrt(constant * alpha[0]) * sigma2;
    }

    public static double[] getConditionedMeanSd(double mu1, double s1, double mu2, double s2,
                                                 double rho, double x1) {
        double mean = mu2 + (s2/s1) * rho * (x1 - mu1);
        double sd = sqrt( (1 - rho * rho) * s2 * s2 );
        return new double[]{mean, sd};
    }

    @Override
    public SortedMap<String, Value> getParams() {
        SortedMap<String, Value> map = new TreeMap<>();
        map.put(muParamName, mu);
        map.put(sigmaParamName, sigma);
        map.put(DRIFT_PARAM, drift);
        map.put(DRIFT_CORR_PARAM, driftCorr);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(muParamName)) mu = value;
        else if (paramName.equals(sigmaParamName)) sigma = value;
        else if (paramName.equals(DRIFT_PARAM)) drift = value;
        else if (paramName.equals(DRIFT_CORR_PARAM)) driftCorr = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }
}
