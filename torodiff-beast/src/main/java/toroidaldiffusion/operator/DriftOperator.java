package toroidaldiffusion.operator;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.operator.kernel.BactrianScaleOperator;

public class DriftOperator extends BactrianScaleOperator {

//    final public Input<Function> driftInput = new Input<>("drift", "the two drift terms.");
    final public Input<Function> driftCorrInput = new Input<>("driftCorr", "the correlation of two drift terms, ranged from -1 to 1.");


    @Override
    public void initAndValidate() {
        boolean success = validateAlpha();
        if (!success)
            throw new IllegalArgumentException("Alpha1 * alpha2 must > alpha3 * alpha3 !");
        super.initAndValidate();
    }

    @Override
    public double proposal() {
        double logProb = super.proposal();
        boolean success = validateAlpha();
        if (!success) return Double.NEGATIVE_INFINITY;
        return logProb;
    }

    private boolean validateAlpha() {
        double[] drift = parameterInput.get().getDoubleValues();
        double corr = driftCorrInput.get().getArrayValue();

        return drift[0] * drift[1] > corr * corr;
    }

}
