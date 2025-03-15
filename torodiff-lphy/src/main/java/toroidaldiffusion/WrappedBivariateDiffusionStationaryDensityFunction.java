package toroidaldiffusion;

import org.apache.commons.math3.analysis.MultivariateFunction;

public class WrappedBivariateDiffusionStationaryDensityFunction implements MultivariateFunction {

    WrappedBivariateDiffusion diff;

    WrappedBivariateDiffusionStationaryDensityFunction(WrappedBivariateDiffusion diff) {
        this.diff = diff;
    }

    @Override
    public double value(double[] doubles) {
        return diff.loglikwndstat(doubles[0], doubles[1]);
    }
}
