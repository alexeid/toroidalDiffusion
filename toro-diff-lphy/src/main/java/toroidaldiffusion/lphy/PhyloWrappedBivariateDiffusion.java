package toroidaldiffusion.lphy;

import lphy.evolution.tree.TimeTree;
import lphy.evolution.tree.TimeTreeNode;
import lphy.graphicalModel.GenerativeDistribution;
import lphy.graphicalModel.ParameterInfo;
import lphy.graphicalModel.RandomVariable;
import lphy.graphicalModel.Value;
import lphy.graphicalModel.types.DoubleArray2DValue;
import lphy.util.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * @author Alexei Drummond
 */
public class PhyloWrappedBivariateDiffusion implements GenerativeDistribution<Map<String, Double[][]>> {

    boolean anglesInRadians = true;

    // ANGLES IN RADIANS FOR THIS IMPLEMENTATION
    double MAX_ANGLE_VALUE = Math.PI * 2.0;

    Value<TimeTree> tree;
    Value<Double[]> mu;
    Value<Double[]> sigma;
    Value<Double[]> alpha;
    Value<Double[][]> y;
    RandomGenerator random;

    public static final String treeParamName = "tree";
    public static final String muParamName = "mu";
    public static final String sigmaParamName = "sigma";
    public static final String alphaParamName = "alpha";
    public static final String y0RateParam = "y0";

    public PhyloWrappedBivariateDiffusion(@ParameterInfo(name = treeParamName, description = "the time tree.") Value<TimeTree> tree,
                                          @ParameterInfo(name = muParamName, description = "the mean of the stationary distribution.") Value<Double[]> mu,
                                          @ParameterInfo(name = sigmaParamName, description = "the two variance terms.") Value<Double[]> sigma,
                                          @ParameterInfo(name = alphaParamName, description = "the three drift terms.") Value<Double[]> alpha,
                                          @ParameterInfo(name = y0RateParam, description = "the value of [phi,psi] angle pairs for each carbon backbone bond of the molecule at the root of the phylogeny.") Value<Double[][]> y) {
        this.tree = tree;
        this.mu = mu;
        this.sigma = sigma;
        this.alpha = alpha;
        this.y = y;
        this.random = RandomUtils.getRandom();
    }

    @Override
    public SortedMap<String, Value> getParams() {
        SortedMap<String, Value> map = new TreeMap<>();
        map.put(treeParamName, tree);
        map.put(muParamName, mu);
        map.put(sigmaParamName, sigma);
        map.put(alphaParamName, alpha);
        map.put(y0RateParam, y);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(treeParamName)) tree = value;
        else if (paramName.equals(muParamName)) mu = value;
        else if (paramName.equals(sigmaParamName)) sigma = value;
        else if (paramName.equals(alphaParamName)) alpha = value;
        else if (paramName.equals(y0RateParam)) y = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

    public RandomVariable<Map<String, Double[][]>> sample() {

        SortedMap<String, Integer> idMap = new TreeMap<>();
        fillIdMap(tree.value().getRoot(), idMap);

        Map<String, Double[][]> tipValues = new StringDoubleArray2DMap();

        WrappedBivariateDiffusion wrappedBivariateDiffusion = new WrappedBivariateDiffusion();

        wrappedBivariateDiffusion.setParameters(mu.value(), alpha.value(), sigma.value());

        traverseTree(tree.value().getRoot(), y, tipValues, wrappedBivariateDiffusion, idMap);

        return new RandomVariable<>("x", tipValues, this);
    }

    private void fillIdMap(TimeTreeNode node, SortedMap<String, Integer> idMap) {
        if (node.isLeaf()) {
            Integer i = idMap.get(node.getId());
            if (i == null) {
                int nextValue = 0;
                for (Integer j : idMap.values()) {
                    if (j >= nextValue) nextValue = j + 1;
                }
                idMap.put(node.getId(), nextValue);
            }
        } else {
            for (TimeTreeNode child : node.getChildren()) {
                fillIdMap(child, idMap);
            }
        }
    }

    private void traverseTree(TimeTreeNode node, Value<Double[][]> nodeState, Map<String, Double[][]> tipValues, WrappedBivariateDiffusion diffusion, Map<String, Integer> idMap) {
        if (node.isLeaf()) {
            tipValues.put(node.getId(), nodeState.value());
        } else {
            for (TimeTreeNode child : node.getChildren()) {

                double branchLength = node.getAge() - child.getAge();

                Double[][] newValue = getNewValue(nodeState.value(), diffusion, branchLength);

                DoubleArray2DValue ns = new DoubleArray2DValue(null, newValue);

                traverseTree(child, ns, tipValues, diffusion, idMap);
            }
        }
    }

    Double[][] getNewValue(Double[][] oldValue, WrappedBivariateDiffusion diffusion, double branchLength) {

        Double[][] newValues = new Double[oldValue.length][oldValue[0].length];

        diffusion.setParameters(branchLength);
        for (int i = 0; i < oldValue.length; i += 1) {

            double[][] samples = diffusion.sampleByRejection(oldValue[i][0], oldValue[i][1], 1);

            newValues[i][0] = samples[0][0];
            newValues[i][1] = samples[0][1];
        }
        return newValues;
    }

}