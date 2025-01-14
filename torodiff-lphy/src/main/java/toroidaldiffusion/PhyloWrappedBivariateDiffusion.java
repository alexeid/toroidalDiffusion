package toroidaldiffusion;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import static java.lang.Math.sqrt;

/**
 * @author Alexei Drummond
 */
public class PhyloWrappedBivariateDiffusion implements GenerativeDistribution<TaxaCharacterMatrix> {

    boolean anglesInRadians = true;

    // ANGLES IN RADIANS FOR THIS IMPLEMENTATION
    double MAX_ANGLE_VALUE = Math.PI * 2.0;

    Value<TimeTree> tree;
    Value<Double[]> mu;
    Value<Double[]> sigma;
    Value<Double[]> drift;
    Value<Number> driftCorr;
    Value<Double[][]> y0;
    RandomGenerator random;
//    Value<Double[]> alpha;

    public PhyloWrappedBivariateDiffusion(@ParameterInfo(name = WrappedNormalConst.treeParamName, description = "the time tree.") Value<TimeTree> tree,
                                          @ParameterInfo(name = WrappedNormalConst.muParamName, description = "the mean of the stationary distribution.") Value<Double[]> mu,
                                          @ParameterInfo(name = WrappedNormalConst.sigmaParamName, description = "the two variance terms.") Value<Double[]> sigma,
                                          @ParameterInfo(name = WrappedNormalConst.DRIFT_PARAM, description = "the two drift terms.") Value<Double[]> drift,
                                          @ParameterInfo(name = WrappedNormalConst.DRIFT_CORR_PARAM, description = "the correlation of two drift terms, ranged from -1 to 1.") Value<Number> driftCorr,
                                          @ParameterInfo(name = WrappedNormalConst.y0RateParam, description = "the value of [phi,psi] angle pairs for each carbon backbone bond of the molecule at the root of the phylogeny.") Value<Double[][]> y0) {
        this.tree = tree;
        this.mu = mu;
        this.sigma = sigma;
        this.drift = drift;
        this.driftCorr = driftCorr;
        this.y0 = y0;
        this.random = RandomUtils.getRandom();

    }

//    @Deprecated
//    public PhyloWrappedBivariateDiffusion(@ParameterInfo(name = treeParamName, description = "the time tree.") Value<TimeTree> tree,
//                                          @ParameterInfo(name = muParamName, description = "the mean of the stationary distribution.") Value<Double[]> mu,
//                                          @ParameterInfo(name = sigmaParamName, description = "the two variance terms.") Value<Double[]> sigma,
//                                          @ParameterInfo(name = alphaParamName, description = "the three drift terms.") Value<Double[]> alpha,
//                                          @ParameterInfo(name = y0RateParam, description = "the value of [phi,psi] angle pairs for each carbon backbone bond of the molecule at the root of the phylogeny.") Value<Double[][]> y) {
//        this.tree = tree;
//        this.mu = mu;
//        this.sigma = sigma;
//        this.alpha = alpha;
//        this.y = y;
//        this.random = RandomUtils.getRandom();
//    }

    @Override
    public SortedMap<String, Value> getParams() {
        SortedMap<String, Value> map = new TreeMap<>();
        map.put(WrappedNormalConst.treeParamName, tree);
        map.put(WrappedNormalConst.muParamName, mu);
        map.put(WrappedNormalConst.sigmaParamName, sigma);
//        map.put(alphaParamName, alpha); // overload
        map.put(WrappedNormalConst.DRIFT_PARAM, drift);
        map.put(WrappedNormalConst.DRIFT_CORR_PARAM, driftCorr);
        map.put(WrappedNormalConst.y0RateParam, y0);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(WrappedNormalConst.treeParamName)) tree = value;
        else if (paramName.equals(WrappedNormalConst.muParamName)) mu = value;
        else if (paramName.equals(WrappedNormalConst.sigmaParamName)) sigma = value;
//        else if (paramName.equals(alphaParamName)) alpha = value; // overload
        else if (paramName.equals(WrappedNormalConst.DRIFT_PARAM)) drift = value;
        else if (paramName.equals(WrappedNormalConst.DRIFT_CORR_PARAM)) driftCorr = value;
        else if (paramName.equals(WrappedNormalConst.y0RateParam)) y0 = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

    @GeneratorInfo(name = "PhyloWrappedBivariateDiffusion", verbClause = "is assumed to have evolved under",
            narrativeName = "phylogenetic bivariate wrapped normal distribution",
            category = GeneratorCategory.PHYLO_LIKELIHOOD, examples = {"simplePhyloWrappedBivariateDiffusion.lphy"},
            description = "The phylogenetic A bivariate wrapped normal distribution distribution.")

    public RandomVariable<TaxaCharacterMatrix> sample() {
        TimeTree timeTree = tree.value();
//        SortedMap<String, Integer> idMap = new TreeMap<>();
//        fillIdMap(timeTree.getRoot(), idMap);
//        Taxa taxa = Taxa.createTaxa(idMap); // TODO Alexei: why not tree.value().getTaxa() ?
        Taxa taxa = timeTree.getTaxa();

        // TODO if y not specified then simulate from equilibrium distribution
        // else just set y0 to given value
        Double[][] y0mat = y0.value();

        int nchar = y0mat.length; // rows
        int length = y0mat[0].length; // cols should be 2, because each element in this array is an angle pair

        if (length != 2)
            throw new RuntimeException("Dimensions of y0 should be L by 2, but found: " + nchar + " by " + length);

        // check if all internal nodes have their IDs
        List<TimeTreeNode> internalNodes = timeTree.getInternalNodes();
        boolean addIntNodeSeq = internalNodes.stream().allMatch(node -> node.getId() != null);
        if (!addIntNodeSeq)
            internalNodes = null; // if internal nodes have no ID, then only tips seqs

        TaxaCharacterMatrix nodeValues = new DihedralAngleAlignment(taxa, nchar, internalNodes);

        WrappedBivariateDiffusion wrappedBivariateDiffusion = new WrappedBivariateDiffusion();

        Double[] alpha_true = getAlphaarr(drift, driftCorr); // should be 3 numbers

        if (alpha_true[0] * alpha_true[1] <= alpha_true[2] * alpha_true[2]) //corr changed to a3
            throw new IllegalArgumentException("Alpha1 * alpha2 must > alpha3 * alpha3 ! But alpha = {" +
                    alpha_true[0] + ", " + alpha_true[1] + ", " + alpha_true[2] +  "} is invalid.");

        /** how alpha is used:
         * double quo = Math.sqrt(sigma.get(0, 0) / sigma.get(1, 0));
         * A.set(0, 0, alpha.get(0, 0));
         * A.set(1, 1, alpha.get(1, 0));
         * A.set(0, 1, alpha.get(2, 0) * quo);
         * A.set(1, 0, alpha.get(2, 0) / quo);
         */

//        Double[] twoDrifts = drift.value();
//        Double[] alpha = new Double[]{twoDrifts[0], twoDrifts[1], driftCorr.value().doubleValue()}; //change to driftCorr -> a3
        wrappedBivariateDiffusion.setParameters(mu.value(), alpha_true, sigma.value());

        // if any internal node id is not null and not empty string,
        // then add its sequence to the alignment.
        traverseTree(timeTree.getRoot(), y0mat, nodeValues, wrappedBivariateDiffusion);

        double[] range = DihedralAngleAlignment.getAngleRange((DihedralAngleAlignment) nodeValues);

        LoggerUtils.log.info("Phi range = [" + range[0] + ", " + range[1] +
                "], Psi range = [" + range[2] + ", " + range[3] + "].");

        return new RandomVariable<>("x", nodeValues, this);
    }

    // compute A given two drifts and their correlation

    /**
     * Reparametrisation: alpha3^2 < alpha2 * alpha3, alpha3 < sqrt(alpha2 * alpha3)
     * @param drift
     * @param driftCorr
     * @return
     */

    public static Double[] getAlphaarr(Value<Double[]> drift, Value<Number> driftCorr) {
        Double[] twoDrifts = drift.value();
        double corr = driftCorr.value().doubleValue();

        double alpha1 = twoDrifts[0].doubleValue();
        double alpha2 = twoDrifts[1].doubleValue();
        double alpha3 = sqrt(alpha1*alpha2) * corr;

        return new Double[]{twoDrifts[0], twoDrifts[1], alpha3};
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

    private void traverseTree(TimeTreeNode node, Double[][] nodeStateArray, TaxaCharacterMatrix<Pair> nodeValues,
                              WrappedBivariateDiffusion diffusion) {

        if (node.isLeaf() || (node.getId() != null && !node.getId().trim().isEmpty())) {
            int nodeIndex = alignmentRowindex(node, nodeValues);

            for (int i = 0; i < nodeStateArray.length; i++) {
                Pair pair = new Pair(nodeStateArray[i][0], nodeStateArray[i][1]);
                // internal nodes
                nodeValues.setState(nodeIndex, i, pair);
            }
        }

        // recursion
        for (TimeTreeNode child : node.getChildren()) {
            double branchLength = node.getAge() - child.getAge();
            Double[][] newValue = getNewValue(nodeStateArray, diffusion, branchLength);
            traverseTree(child, newValue, nodeValues, diffusion);
        }
    }

    private int alignmentRowindex(TimeTreeNode node, TaxaCharacterMatrix<Pair> alignment) {
        if (node.isLeaf()) {
            String id = node.getId();
            int index = alignment.getTaxa().indexOfTaxon(id);
            return index;
        }
        return node.getIndex();
    }

    Double[][] getNewValue(Double[][] oldValue, WrappedBivariateDiffusion diffusion, double branchLength) {

        Double[][] newValues = new Double[oldValue.length][oldValue[0].length];

        diffusion.setParameters(branchLength);
        for (int i = 0; i < oldValue.length; i += 1) {

            try {
                double[][] samples = diffusion.sampleByRejection(oldValue[i][0], oldValue[i][1], 1);

                newValues[i][0] = samples[0][0];
                newValues[i][1] = samples[0][1];
            } catch (MathIllegalStateException e) {
                throw new IllegalStateException("Cannot calculate new value for " + oldValue[i][0] +
                        " and " + oldValue[i][1] + ", given branch length = " + branchLength +
                        "\nmu = [" + diffusion.mu.get(0) + ", " + diffusion.mu.get(1) + "]" +
                        "\nsigma = [" + diffusion.sigma.get(0) + ", " + diffusion.sigma.get(1) + "]" +
                        "\nalpha = [" + diffusion.alpha.get(0) + ", " + diffusion.alpha.get(1) + ", " + diffusion.alpha.get(2) + "]");
            }
        }

        return newValues;
    }

    public Value<TimeTree> getTree() {
        return tree;
    }

    public Value<Double[][]> getY0() {
        return y0;
    }

}