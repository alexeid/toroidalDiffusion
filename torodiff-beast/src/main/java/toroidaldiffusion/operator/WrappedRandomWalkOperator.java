package toroidaldiffusion.operator;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

import java.text.DecimalFormat;

public class WrappedRandomWalkOperator extends Operator {

    final public Input<Double> windowSizeInput =
            new Input<>("windowSize", "the size of the window both up and down when using uniform interval OR standard deviation when using Gaussian", Input.Validate.REQUIRED);
    final public Input<RealParameter> parameterInput =
            new Input<>("parameter", "the parameter to operate a random walk on.", Input.Validate.REQUIRED);
    final public Input<Boolean> useGaussianInput =
            new Input<>("useGaussian", "Use Gaussian to move instead of uniform interval. Default false.", false);

    public static final double ANGLE_LOWER = 0.0;
    public static final double ANGLE_UPPER = 2 * Math.PI;

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

        // Wrap the new value to be within 0 to 2π
        newValue = wrapAngles(newValue);

        // If the new value is the same as the current value, reject the proposal
        if (newValue == value) {
            return Double.NEGATIVE_INFINITY;
        }

        // Set the new value
        param.setValue(i, newValue);

        return 0.0; // Hastings ratio is 0.0 (log(1))
    }

    // TODO range [0, 2pi)
    //  > 2 * Math.PI or >=
    public static double wrapAngles(double angle) {
        while (angle >= ANGLE_UPPER) {
            angle -= 2 * Math.PI;
        }
        while (angle < ANGLE_LOWER) {
            angle += 2 * Math.PI;
        }
        return angle;
    }

    public double getCoercableParameterValue() {
        return windowSize;
    }

    public void setCoercableParameterValue(double v) {
        windowSize = v;
    }

    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double delta = calcDelta(logAlpha);

        delta += Math.log(windowSize);
        windowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = windowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }

}
