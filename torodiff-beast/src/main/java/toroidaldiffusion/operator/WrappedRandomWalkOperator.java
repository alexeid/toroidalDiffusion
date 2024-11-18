package toroidaldiffusion.operator;

import beast.base.inference.operator.RealRandomWalkOperator;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

public class WrappedRandomWalkOperator extends RealRandomWalkOperator {

    private double windowSize = 1;
    private boolean useGaussian;

    public void initAndValidate() {
        windowSize = windowSizeInput.get();
        useGaussian = useGaussianInput.get();
    }

    public double proposal() {
        RealParameter param = (RealParameter) InputUtil.get(parameterInput, this);

        int i = Randomizer.nextInt(param.getDimension());
        double value = param.getValue(i);
        double newValue;

        // Perturb the value
        if (useGaussian) {
            newValue = value + Randomizer.nextGaussian() * windowSize;
        } else {
            newValue = value + Randomizer.nextDouble() * 2 * windowSize - windowSize;
        }

        // Wrap the new value to be within -π to +π
        newValue = wrapAngles(newValue);

        // If the new value is the same as the current value, reject the proposal
        if (newValue == value) {
            return Double.NEGATIVE_INFINITY;
        }

        // Set the new value
        param.setValue(i, newValue);

        return 0.0; // Hastings ratio is 0.0 (log(1))

        /* compiled code */
    }

    // TODO range [0, 2pi)
    //  > 2 * Math.PI or >=
    public static double wrapAngles(double angle) {
        while (angle >= 2 * Math.PI) {
            angle -= 2 * Math.PI;
        }
        while (angle < 0) {
            angle += 2 * Math.PI;
        }
        return angle;
    }

    public double getCoercableParameterValue() {
        /* compiled code */
        return windowSize;
    }

    public void setCoercableParameterValue(double v) {
        windowSize = v;
        /* compiled code */ }

    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double delta = calcDelta(logAlpha);

        delta += Math.log(windowSize);
        windowSize = Math.exp(delta);
        /* compiled code */ }

}
