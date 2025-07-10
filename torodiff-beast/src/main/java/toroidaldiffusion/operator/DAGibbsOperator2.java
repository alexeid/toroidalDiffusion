package toroidaldiffusion.operator;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import toroidaldiffusion.WrappedBivariateDiffusion;
import toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;

public class DAGibbsOperator2 extends Operator {

    final public Input<PhyloWrappedBivariateDiffusion> PhyloWrappedBivariateDiffusionInput = new Input<>("phyloWrappedBivariateDiffusion",
            "Likelihood calculation class",
            Input.Validate.REQUIRED);
    final public Input<RealParameter> parameterInput =
            new Input<>("parameter", "the parameter to operate a random walk on.",
                    Input.Validate.REQUIRED);
    final public Input<Double> varianceInflationInput =
            new Input<>("varianceInflation", "", 1.0);

    private DihedralAngleTreeModel daTreeModel;
    private WrappedBivariateDiffusion diff;
    private DihedralAngleGibbsSampler2 sampler;
    private double varianceInflation = 1.0;

    @Override
    public void initAndValidate() {
        this.daTreeModel = (DihedralAngleTreeModel) PhyloWrappedBivariateDiffusionInput.get().getDaTreeModel();
        this.diff = PhyloWrappedBivariateDiffusionInput.get().getDiff();
        this.sampler = new DihedralAngleGibbsSampler2();//need to set parameters
        this.varianceInflation = varianceInflationInput.get();

        sampler.setDiff(diff);

//        super.initAndValidate();
    }

    @Override
    public double proposal() {
        // 1. Select a random internal node
        //+ getInternalNodeCounts
        TreeInterface tree = daTreeModel.getTree();

        final int nodeNr = tree.getLeafNodeCount() + Randomizer.nextInt(tree.getInternalNodeCount());
        final Node node = daTreeModel.getNode(nodeNr);

        RealParameter internalNodesAngles = parameterInput.get();

        // 2. calculate the index of this internal node (subtract leaf node count)
        int paramIndex = nodeNr - daTreeModel.getLeafNodeCount();

        // Calculate start and end indices in the parameter array
        // Assuming each internal node has siteCount * 2 values (phi and psi for each site)
        int siteCount = daTreeModel.getSiteCount();
        int startInclusive = paramIndex * siteCount * 2;
        int endExclusive = startInclusive + siteCount * 2;

        //Validation: endIndex > internalNodesAngles.length then out of boundary (not ==, because endIndex exclusive index = last index + 1
        if (endExclusive > internalNodesAngles.getDimension()) {
            throw new IllegalArgumentException("Parameter array index out of bounds: calculated end index "
                    + endExclusive + " exceeds parameter dimension "
                    + internalNodesAngles.getDimension());
        }

        // 4. Perform Gibbs sampling to get new angles
        double[] proposedAngles = sampler.gibbsSampling(node, daTreeModel, varianceInflation);

//        // Accept the proposal
        for (int i = 0; i < proposedAngles.length; i++) {
            internalNodesAngles.setValue(startInclusive + i, proposedAngles[i]);
        }

        // 5. Calculate Hastings ratio
//        double logHastingsRatio = sampler.calculateLogHastingsRatio(currentAngles, proposedAngles, daTreeModel);
        double logHastingsRatio = sampler.getLogHastingsratio();

        return logHastingsRatio;

    }
}