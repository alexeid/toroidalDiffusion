package toroidaldiffusion;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.SortedMap;
import java.util.TreeMap;

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
    Value<Double[]> alpha;
    Value<Double[][]> y;
    RandomGenerator random;

    public static final String treeParamName = "tree";
    public static final String muParamName = "mu";
    public static final String sigmaParamName = "sigma";
    public static final String alphaParamName = "alpha";
    public static final String y0RateParam = "y0";
    public static final String ADD_INT_NODE_SEQ = "addIntNodeSeq";

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

    @GeneratorInfo(name = "PhyloWrappedBivariateDiffusion", verbClause = "is assumed to have evolved under",
            narrativeName = "phylogenetic bivariate wrapped normal distribution",
            category = GeneratorCategory.PHYLO_LIKELIHOOD, examples = {"simplePhyloWrappedBivariateDiffusion.lphy"},
            description = "The phylogenetic A bivariate wrapped normal distribution distribution.")
    public RandomVariable<TaxaCharacterMatrix> sample() {
        TimeTree timeTree = tree.value();
        SortedMap<String, Integer> idMap = new TreeMap<>();
        fillIdMap(timeTree.getRoot(), idMap);
        Taxa taxa = Taxa.createTaxa(idMap); // TODO Alexei: why not tree.value().getTaxa() ?

        Double[][] y0 = y.value();

        int nchar = y0.length;
        int length = y0[0].length; // should be 2, because each element in this array is an angle pair

        if (length != 2)
            throw new RuntimeException("Dimensions of y0 should be L by 2, but found: " + nchar + " by " + length);

        // check if all internal nodes have their IDs
        boolean addIntNodeSeq = timeTree.getInternalNodes().stream().allMatch(node -> node.getId() != null);
        TaxaCharacterMatrix nodeValues = new DihedralAngleAlignment(taxa, nchar, addIntNodeSeq);

        WrappedBivariateDiffusion wrappedBivariateDiffusion = new WrappedBivariateDiffusion();

        wrappedBivariateDiffusion.setParameters(mu.value(), alpha.value(), sigma.value());

        // if any internal node id is not null and not empty string,
        // then add its sequence to the alignment.
        traverseTree(timeTree.getRoot(), y0, nodeValues, wrappedBivariateDiffusion);

        return new RandomVariable<>("x", nodeValues, this);
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

        if (node.isLeaf()) {
            Taxa taxa = nodeValues.getTaxa();
            int taxonIndex = taxa.indexOfTaxon(node.getId());

            for (int i = 0; i < nodeStateArray.length; i++) {
                Pair pair = new Pair(nodeStateArray[i][0], nodeStateArray[i][1]);
                // tips
                nodeValues.setState(taxonIndex, i, pair);
            }

        } else {
            for (TimeTreeNode child : node.getChildren()) {

                double branchLength = node.getAge() - child.getAge();

                Double[][] newValue = getNewValue(nodeStateArray, diffusion, branchLength);
                // if id not null and not empty string, then add internal nodes sequences to the alignment.
                if (node.getId() != null && !node.getId().trim().isEmpty()) {
                    int nodeIndex = node.getIndex();
                    for (int i = 0; i < newValue.length; i++) {
                        Pair pair = new Pair(newValue[i][0], newValue[i][1]);
                        // internal nodes
                        nodeValues.setState(nodeIndex, i, pair);
                    }
                }

                traverseTree(child, newValue, nodeValues, diffusion);
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