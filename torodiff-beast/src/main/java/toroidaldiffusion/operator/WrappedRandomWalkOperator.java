package toroidaldiffusion.operator;

import beast.base.core.Input;
import beast.base.inference.operator.RealRandomWalkOperator;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

public class WrappedRandomWalkOperator extends RealRandomWalkOperator {
    final public Input<Double> windowSizeInput =
            new Input<>("windowSize", "the size of the window both up and down when using uniform interval OR standard deviation when using Gaussian", Input.Validate.REQUIRED);
    final public Input<RealParameter> parameterInput =
            new Input<>("parameter", "the parameter to operate a random walk on.", Input.Validate.REQUIRED);
    final public Input<Boolean> useGaussianInput =
            new Input<>("useGaussian", "Use Gaussian to move instead of uniform interval. Default false.", false);

    double windowSize;
    boolean useGaussian;

    public WrappedRandomWalkOperator() {
    }

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

    private double wrapAngles(double angle) {
        while (angle > Math.PI) {
            angle -= 2 * Math.PI;
        }
        while (angle < -Math.PI) {
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
