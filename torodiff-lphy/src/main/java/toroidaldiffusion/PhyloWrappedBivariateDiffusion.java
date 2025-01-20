package toroidaldiffusion;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.random.RandomGenerator;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import static java.lang.Math.sqrt;
import static toroidaldiffusion.WrappedNormalConst.*;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PhyloWrappedBivariateDiffusion implements GenerativeDistribution<TaxaCharacterMatrix> {

    boolean anglesInRadians = true;

    Value<TimeTree> tree;
    Value<Number[]> mu;
    Value<Number[]> sigma;
    Value<Number[]> drift;
    Value<Number> driftCorr;
    Value<Integer> l;
//    Value<Double[][]> y0;

    RandomGenerator random;
//    Value<Double[]> alpha;

    public PhyloWrappedBivariateDiffusion(@ParameterInfo(name = treeParamName, description = "the time tree.") Value<TimeTree> tree,
                                          @ParameterInfo(name = muParamName, description = "the mean of the stationary distribution.") Value<Number[]> mu,
                                          @ParameterInfo(name = sigmaParamName, description = "the two variance terms.") Value<Number[]> sigma,
                                          @ParameterInfo(name = DRIFT_PARAM, description = "the two drift terms.") Value<Number[]> drift,
                                          @ParameterInfo(name = DRIFT_CORR_PARAM, description = "the correlation of two drift terms, ranged from -1 to 1.") Value<Number> driftCorr,
                                          @ParameterInfo(name = LParamName, description = "the number of pairs of angles (sites) to be simulated.") Value<Integer> l) {
//                                          @ParameterInfo(name = WrappedNormalConst.y0RateParam, description = "the value of [phi,psi] angle pairs for each carbon backbone bond of the molecule at the root of the phylogeny.") Value<Double[][]> y0) {
        this.tree = tree;
        this.mu = mu;
        this.sigma = sigma;
        this.drift = drift;
        this.driftCorr = driftCorr;
        this.l = l;
//        this.y0 = y0;
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
        map.put(treeParamName, tree);
        map.put(muParamName, mu);
        map.put(sigmaParamName, sigma);
//        map.put(alphaParamName, alpha); // overload
        map.put(DRIFT_PARAM, drift);
        map.put(DRIFT_CORR_PARAM, driftCorr);
        map.put(LParamName, l);
//        map.put(y0RateParam, y0);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(treeParamName)) tree = value;
        else if (paramName.equals(muParamName)) mu = value;
        else if (paramName.equals(sigmaParamName)) sigma = value;
//        else if (paramName.equals(alphaParamName)) alpha = value; // overload
        else if (paramName.equals(DRIFT_PARAM)) drift = value;
        else if (paramName.equals(DRIFT_CORR_PARAM)) driftCorr = value;
        else if (paramName.equals(LParamName)) l = value;
//        else if (paramName.equals(y0RateParam)) y0 = value;
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

        int nchar = l.value(); // sites

        // check if all internal nodes have their IDs
        List<TimeTreeNode> internalNodes = timeTree.getInternalNodes();
        boolean addIntNodeSeq = internalNodes.stream().allMatch(node -> node.getId() != null);
        if (!addIntNodeSeq)
            internalNodes = null; // if internal nodes have no ID, then only tips seqs

        TaxaCharacterMatrix nodeValues = new DihedralAngleAlignment(taxa, nchar, internalNodes);

        WrappedBivariateDiffusion wrappedBivariateDiffusion = new WrappedBivariateDiffusion();

        double[] alpha_true = getAlphaArr(drift, driftCorr); // should be 3 numbers

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
        double[] muArr = ValueUtils.doubleArrayValue(mu);
        wrappedBivariateDiffusion.setParameters(muArr, alpha_true, ValueUtils.doubleArrayValue(sigma));

        // root sequences y0 should be simulated from equilibrium distribution
        // sampling method both on root sequences and other internal nodes should be same,
        // method 1: rejection sampling at mu
        // TODO wrong, sampling from stationary dist does not require time, but sampleByRejection does
//        Double[][] y0 = simulateRootSeqs(wrappedBivariateDiffusion, muArr[0], muArr[1], nchar);
        // method 2: sampling from WN(mu, 1/2 A^-1 Sigma)
        Double[][] y0 = simulateRootSeqs2();

        // if any internal node id is not null and not empty string,
        // then add its sequence to the alignment.
        traverseTree(timeTree.getRoot(), y0, nodeValues, wrappedBivariateDiffusion);

        double[] range = DihedralAngleAlignment.getAngleRange((DihedralAngleAlignment) nodeValues);

        LoggerUtils.log.info("Phi range = [" + range[0] + ", " + range[1] +
                "], Psi range = [" + range[2] + ", " + range[3] + "].");

        return new RandomVariable<>("x", nodeValues, this);
    }

    /**
     * rejection sampling at mu
     * @param wrappedBiDif
     * @param muPhi
     * @param muPsi
     * @param nsamples
     * @return
     */
    public static Double[][] simulateRootSeqs(WrappedBivariateDiffusion wrappedBiDif,
                                       double muPhi, double muPsi, int nsamples){
        //TODO this is wrong, sampling from stationary dist does not require time, but sampleByRejection does
        double[][] y0 = wrappedBiDif.sampleByRejection(muPhi, muPsi, nsamples);
        return Arrays.stream(y0).map(
                row -> Arrays.stream(row).boxed().toArray(Double[]::new)
        ).toArray(Double[][]::new);
    }

    /**
     * simulate root sequences y0 from equilibrium distribution WN(mu, 1/2 A^-1 Sigma)
     * @return y0
     */
    public Double[][] simulateRootSeqs2(){
        int nchar = l.value(); // number of sites
        Double[][] y0 = new Double[nchar][2];
        WrappedBivariateNormal equilDist = new WrappedBivariateNormal(mu, sigma, drift, driftCorr);
        for (int i = 0; i < nchar; i++) {
            Double[] pair = equilDist.sample().value();
            y0[i] = pair;
        }
        return y0;
    }

    // compute A given two drifts and their correlation

    /**
     * Reparametrisation: alpha3^2 < alpha1 * alpha2, alpha3 < sqrt(alpha1 * alpha2)
     * @param drift
     * @param driftCorr
     * @return
     */
    public static double[] getAlphaArr(Value<Number[]> drift, Value<Number> driftCorr) {
        double[] twoDrifts = ValueUtils.doubleArrayValue(drift);
        double corr = ValueUtils.doubleValue(driftCorr);

        double alpha1 = twoDrifts[0];
        double alpha2 = twoDrifts[1];
        double alpha3 = sqrt(alpha1*alpha2) * corr;

        return new double[]{alpha1, alpha2, alpha3};
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
                throw new IllegalStateException("Cannot calculate new (phi, psi) given old pair (" + oldValue[i][0] +
                        ", " + oldValue[i][1] + "), where branch length = " + branchLength +
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

//    public Value<Double[][]> getY0() {
//        return y0;
//    }

    public static void main(String[] args) throws IOException {

        final int NSamples = 100000;

        WrappedBivariateDiffusion diff = new WrappedBivariateDiffusion();
        double[] muarr = {Math.PI * 0.65, Math.PI * 0.8}; // mean of the diffusion
        double[] sigmaarr = {1.5, 1.5}; // variance term
        double[] alphaarr = {1.0, 0.5, 0.5}; // drift term
        // before this method time is set to 1.0 as default
        diff.setParameters(muarr, alphaarr, sigmaarr); // set the diffusion parameters

        System.out.println("Stationary dist logP = " + diff.loglikwndstat(muarr[0], muarr[1])); // calculate the stationary density of the mean
        //1. rejection sampling at mu
        Double[][] y0 = simulateRootSeqs(diff, muarr[0], muarr[1], NSamples);

        String filename = "sampleByRejectionAtMu.txt";
        PrintWriter writer = new PrintWriter(new FileWriter(filename));
        writer.println("phi\tpsi\tlogP\tdensity");
        for (int i = 0; i < y0.length; i++) {
            double phi = y0[i][0];
            double psi = y0[i][1];

            double logP = diff.loglikwndstat(phi, psi);
            writer.println(phi + "\t" + psi + "\t" + logP + "\t" + Math.exp(logP)); // calculate the transition density of the point (0.0, 0.0) transitioning to (1.0, 1.0) in time t=1.0
        }
        writer.flush();
        writer.close();

        // 2. sampling from WN(mu, 1/2 A^-1 Sigma)
//        Double[][] y0 = simulateRootSeqs2();


        // prove t does not affect the loglikwndstat code
//        String filename = "testStatDistJava.txt";
//        PrintWriter writer = new PrintWriter(new FileWriter(filename));
//        writer.println("phi\tpsi\tlogP\tdensity");
//        for (int i = 0; i < 1000; i++) {
//            double phi = muarr[0];
//            double psi = muarr[1];
//
//            diff.setParameters(i * 0.5);
//            double logP = diff.loglikwndstat(phi, psi);
//            writer.println(phi + "\t" + psi + "\t" + logP + "\t" + Math.exp(logP) + "\t" + i * 0.5); // calculate the transition density of the point (0.0, 0.0) transitioning to (1.0, 1.0) in time t=1.0
//        }
//        writer.flush();
//        writer.close();

    }
}