package toroidaldiffusion;

import lphy.base.distribution.ParametricDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.SortedMap;
import java.util.TreeMap;

import static java.lang.Math.sqrt;
import static toroidaldiffusion.WrappedNormalConst.MAX_ANGLE_VALUE;

public class WrappedBivariateNormal extends ParametricDistribution<Double[]> {

    Value<Number[]> mu;
    Value<Number[]> sigma;
    Value<Number[]> drift;
    Value<Number> driftCorr;

    NormalDistribution nD1;
    NormalDistribution nD2;

    public WrappedBivariateNormal(@ParameterInfo(name = WrappedNormalConst.muParamName, description = "the mean of the stationary distribution.")
                                  Value<Number[]> mu,
                                  @ParameterInfo(name = WrappedNormalConst.sigmaParamName, description = "the two variance terms.")
                                  Value<Number[]> sigma,
                                  @ParameterInfo(name = WrappedNormalConst.DRIFT_PARAM, description = "the two drift terms : alpha1 and alpha2.")
                                  Value<Number[]> drift,
                                  @ParameterInfo(name = WrappedNormalConst.DRIFT_CORR_PARAM, description = "the correlation of two drift terms, " +
                                          "ranged from -1 to 1 excluding -1 and 1, but not alpha3 where alpha3 = sqrt(alpha1*alpha2) * corr.")
                                  Value<Number> driftCorr) {
        super(); // getRandom
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

    @GeneratorInfo(name = "WrappedBivariateNormal", verbClause = "has", narrativeName = "normal prior",
            category = GeneratorCategory.PRIOR, description = "The wrapped bivariate normal distribution (Golden et al. 2017).")
    @Override
    public RandomVariable<Double[]> sample() {
        double x1 = nD1.sample();

        double mu1 = ValueUtils.doubleArrayValue(mu)[0];
        double mu2 = ValueUtils.doubleArrayValue(mu)[1];
        // transfer alpha1,2,3 into rho
        double[] alpha = WrappedNormalUtils.getAlphaArr(drift, driftCorr);
        double rho = getRho(alpha);
        // transfer sigma1,2 into s1,2 for bivariate normal
        double s1 = getS1(ValueUtils.doubleArrayValue(sigma)[0], alpha);
        double s2 = getS2(ValueUtils.doubleArrayValue(sigma)[1], alpha);

        // sample from the second Normal conditioned on the first Normal
        // X2 <- rnorm(n, mu2 + (s2/s1) * rho * (X1 - mu1), sqrt((1 - rho^2)*s2^2))
        double[] mean_sd = getConditionedMeanSd(mu1, s1, mu2, s2, rho, x1);
        nD2 = new NormalDistribution(random, mean_sd[0], mean_sd[1], NormalDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
        double x2 = nD2.sample();

        double x1w = ToroidalUtils.moduloInR(x1, MAX_ANGLE_VALUE);
        double x2w = ToroidalUtils.moduloInR(x2, MAX_ANGLE_VALUE);

        return new RandomVariable<>("bivariate", new Double[]{x1w, x2w}, this);
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

    // X2 <- rnorm(n, mu2 + (s2/s1) * rho * (X1 - mu1), sqrt((1 - rho^2)*s2^2))
    public static double[] getConditionedMeanSd(double mu1, double s1, double mu2, double s2,
                                                 double rho, double x1) {
        double mean = mu2 + (s2/s1) * rho * (x1 - mu1);
        double sd = sqrt( (1 - rho * rho) * s2 * s2 );
        return new double[]{mean, sd};
    }

    @Override
    public SortedMap<String, Value> getParams() {
        SortedMap<String, Value> map = new TreeMap<>();
        map.put(WrappedNormalConst.muParamName, mu);
        map.put(WrappedNormalConst.sigmaParamName, sigma);
        map.put(WrappedNormalConst.DRIFT_PARAM, drift);
        map.put(WrappedNormalConst.DRIFT_CORR_PARAM, driftCorr);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(WrappedNormalConst.muParamName)) mu = value;
        else if (paramName.equals(WrappedNormalConst.sigmaParamName)) sigma = value;
        else if (paramName.equals(WrappedNormalConst.DRIFT_PARAM)) drift = value;
        else if (paramName.equals(WrappedNormalConst.DRIFT_CORR_PARAM)) driftCorr = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

    public static void main(String[] args) throws IOException {

        final int NSamples = 100000;
        double[] muarr = {Math.PI * 0.65, Math.PI * 0.8}; // mean of the diffusion
        double[] sigmaarr = {1.5, 1.5}; // variance term
        double[] alphaarr = {1.0, 0.5, 0.5}; // drift term
        double driftCorr = WrappedNormalUtils.getDriftCorr(alphaarr[0], alphaarr[1], alphaarr[2]);

        WrappedBivariateNormal wrappedBiNorm = new WrappedBivariateNormal(
                new Value<>("mu", new Double[]{muarr[0], muarr[1]}),
                new Value<>("sigma", new Double[]{sigmaarr[0], sigmaarr[1]}),
                new Value<>("drift", new Double[]{alphaarr[0], alphaarr[1]}),
                new Value<>("driftCorr", driftCorr)
        );

        //2. sampling from WN(mu, 1/2 A^-1 Sigma)

        String filename = "WrappedBivariateNormalJava.txt";
        PrintWriter writer = new PrintWriter(new FileWriter(filename));
        writer.println("phi\tpsi\tlogP\tdensity");
        for (int i = 0; i < NSamples; i++) {
            Double[] pair = wrappedBiNorm.sample().value();

            double phi = pair[0];
            double psi = pair[1];

//            double logP = diff.loglikwndstat(phi, psi);
            writer.println(phi + "\t" + psi); // calculate the transition density of the point (0.0, 0.0) transitioning to (1.0, 1.0) in time t=1.0
        }
        writer.flush();
        writer.close();

    }
}
