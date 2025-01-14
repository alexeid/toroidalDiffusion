package toroidaldiffusion;

import lphy.base.evolution.continuous.PhyloOU;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;

import static toroidaldiffusion.ToroidalUtils.wrapToMaxAngle;
import static toroidaldiffusion.WrappedNormalConst.y0RateParam;

/**
 * Created by adru001 on 2/02/20.
 */
public class PhyloCircularOU extends PhyloOU {

    public PhyloCircularOU(@ParameterInfo(name = treeParamName, description = "the time tree.") Value<TimeTree> tree,
                           @ParameterInfo(name = "variance", description = "the variance of the underlying Brownian process.") Value<Double> variance,
                           @ParameterInfo(name = "theta", description = "the 'optimal' value that the long-term process is centered around.") Value<Double> theta,
                           @ParameterInfo(name = alphaParamName, description = "the drift term that determines the rate of drift towards the optimal value.") Value<Double> alpha,
                           @ParameterInfo(name = y0RateParam, description = "the value of continuous trait at the root.") Value<Double> y0) {

        super(tree, variance, theta, alpha, y0, null);
    }

    @Override
    protected double handleBoundaries(double rawValue) {
        return wrapToMaxAngle(rawValue, 2*Math.PI);
    }
}
