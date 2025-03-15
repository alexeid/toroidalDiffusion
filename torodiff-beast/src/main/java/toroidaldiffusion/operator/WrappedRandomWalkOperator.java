package toroidaldiffusion.operator;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import toroidaldiffusion.ToroidalUtils;

import java.text.DecimalFormat;

public class WrappedRandomWalkOperator extends Operator {

    final public Input<Double> windowSizeInput =
            new Input<>("windowSize", "the size of the window both up and down when using uniform interval OR standard deviation when using Gaussian",
                    Input.Validate.REQUIRED);
    final public Input<RealParameter> parameterInput =
            new Input<>("parameter", "the parameter to operate a random walk on.",
                    Input.Validate.REQUIRED);
    final public Input<Boolean> useGaussianInput =
            new Input<>("useGaussian", "Use Gaussian to move instead of uniform interval. Default false.", false);

//    final public Input<Tree> treeInput =
//            new Input<>("tree", "beast.tree on which this operation is performed",
//                    Input.Validate.REQUIRED);

    public static final double ANGLE_LOWER = 0.0;
    public static final double ANGLE_UPPER = 2 * Math.PI;

    private double windowSize = 1;
    private boolean useGaussian;

//    private Tree tree;

    public void initAndValidate() {
        windowSize = windowSizeInput.get();
        useGaussian = useGaussianInput.get();
//        tree = treeInput.get();
    }

    // need to move the pair of angles each time
    public double proposal() {
        RealParameter param = (RealParameter) InputUtil.get(parameterInput, this);

        int i = Randomizer.nextInt(param.getDimension());
        // (phi, psi) index starts (0, 1), so if i is even then j = i+1, else j = i-1
        int j = i % 2 == 0 ? i + 1 : i - 1;
        double val1 = param.getValue(i);
        double val2 = param.getValue(j);
        double newVal1 = getNewValue(val1);
        double newVal2 = getNewValue(val2);
        // If the new value is the same as the current value, reject the proposal
        if (newVal1 == val1 && newVal2 == val2)
            return Double.NEGATIVE_INFINITY;

        // Set the new pair
        param.setValue(i, newVal1);
        param.setValue(j, newVal2);

//        // assuming nodeIndex is Nr
//        int nodeIndex = (int) Math.floor((double) i / param.getMinorDimension1());
//        int nr = nodeIndex + tree.getLeafNodeCount();
//        if (nr > tree.getNodeCount())
//            throw new IllegalArgumentException("Node index (" + nr + ") must < nodes count (" + tree.getNodeCount() + ") ! ");
//        Node internal = tree.getNode(nr);
//        internal.makeDirty(Tree.IS_DIRTY);

        return 0.0; // Hastings ratio is 0.0 (log(1))
    }

    private double getNewValue(double value) {
        double newValue;
        // Perturb the value
        if (useGaussian) {
            newValue = value + Randomizer.nextGaussian() * windowSize;
        } else {
            newValue = value + Randomizer.nextDouble() * 2 * windowSize - windowSize;
        }

        // Wrap the new value to be within 0 to 2π
        newValue = ToroidalUtils.wrapToMaxAngle(newValue);
        return newValue;
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
