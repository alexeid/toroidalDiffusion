package toroidaldiffusion.operator;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import toroidaldiffusion.WrappedBivariateDiffusion;
import toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;

public class DAGibbsOperator extends Operator {

    final public Input<PhyloWrappedBivariateDiffusion> PhyloWrappedBivariateDiffusionInput = new Input<>("phyloWrappedBivariateDiffusion",
            "Likelihood calculation class",
            Input.Validate.REQUIRED);
    final public Input<RealParameter> parameterInput =
            new Input<>("parameter", "the parameter to operate a random walk on.",
                    Input.Validate.REQUIRED);

    private DihedralAngleTreeModel daTreeModel;
    private WrappedBivariateDiffusion diff;
    private DihedralAngleGibbsSampler sampler;


    @Override
    public void initAndValidate() {
        this.daTreeModel = (DihedralAngleTreeModel) PhyloWrappedBivariateDiffusionInput.get().getDaTreeModel();
        this.diff = PhyloWrappedBivariateDiffusionInput.get().getDiff();
        this.sampler = new DihedralAngleGibbsSampler();//need to set parameters
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

        RealParameter internalNodesAngles = (RealParameter) InputUtil.get(parameterInput, this);

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

        // Store the current angles before proposing new ones
//        double[] currentAngles = new double[siteCount * 2];
//        for (int i = 0; i < currentAngles.length; i++) {
//            currentAngles[i] = internalNodesAngles.getValue(startInclusive + i);
//        }

        // 3. Prepare the sampler with this node
//        sampler.setNode();

        // 4. Perform Gibbs sampling to get new angles
        double[] proposedAngles = sampler.gibbsSampling(node, daTreeModel);

        // 5. Calculate Hastings ratio
//        double logHastingsRatio = sampler.calculateLogHastingsRatio(currentAngles, proposedAngles, daTreeModel);
        double logHastingsRatio = sampler.getLogHastingsratio();

        // 6. Apply Metropolis-Hastings acceptance criteria
        if (logHastingsRatio >= 0.0 || Math.log(Randomizer.nextDouble()) < logHastingsRatio) {
            // Accept the proposal
            for (int i = 0; i < proposedAngles.length; i++) {
                internalNodesAngles.setValue(startInclusive + i, proposedAngles[i]);
            }
            return logHastingsRatio;
        } else {
            // Reject the proposal - no changes needed as we haven't modified the parameter yet
            return Double.NEGATIVE_INFINITY;  // Return -Infinity to signal rejection

        }
        //todo: what to return? since we have calculated ratio for a list of angles
    }
}