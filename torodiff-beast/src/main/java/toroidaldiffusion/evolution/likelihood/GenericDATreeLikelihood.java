package toroidaldiffusion.evolution.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import toroidaldiffusion.evolution.tree.DATreeModel;

import java.util.List;
import java.util.Random;

@Description("Generic data augmentation tree likelihood for data in all nodes " +
        "given a beast // 2 values")
public class GenericDATreeLikelihood extends Distribution {

    final public Input<DATreeModel> daTreeModelInput = new Input<>("daTreeModel",
            "States in all nodes with the beast.tree", Input.Validate.REQUIRED);

    /**
     * {@link DATreeModel}
     */
    protected DATreeModel daTreeModel;
    /**
     * {@link TreeInterface}
     */
    protected TreeInterface tree;
    /**
     * tip values
     */
    protected RealParameter tipValuesParam;
    /**
     * internal nodes values:
     * minordimension / 2 is the number of sites,
     * dimension / minordimension is the node index.
     */
    protected RealParameter internalNodesParam;




    @Override
    public void initAndValidate() {
        // tree model
        daTreeModel = daTreeModelInput.get();
        // tree
        tree = daTreeModel.getTree();
        // data
        tipValuesParam = daTreeModel.getTipsValuesParam();
        internalNodesParam = daTreeModel.getInternalNodesValuesParam();

        if (internalNodesParam.getMinorDimension2() != tree.getInternalNodeCount())
            throw new IllegalArgumentException("The internal node sequences parameter minor dim 2 must = number of internal nodes !");
    }

    /**
     * This assumes internal nodes sequences stored by the order of Nr
     * in the {@link RealParameter} {@link #internalNodesParam}.
     *
     * @param node can be any node
     * @return  whether the internal node sequences mapping to the node given a node index have been changed or not,
     *          where dirty means its value was changed.
     */
    public boolean isInternalNodeSeqDirty(final Node node) {
        // TODO assuming tip sequences are fixed
        if (node.isLeaf()) return false;

        // assuming internalNodeIndex is Nr, and index starts from 1 to the number of nodes.
        int nTips = tree.getLeafNodeCount();
        final int nodeNr = node.getNr();
        // the index to be used for internalNodesParam
        int arrI = nodeNr - nTips;

        // minor dim 2 = number of internal nodes
        if (arrI >= internalNodesParam.getMinorDimension2())
            throw new IllegalArgumentException("Internal node index (" + nodeNr + ") - nTips (" +
                    tree.getLeafNodeCount() + ") cannot > the internal node sequences parameter 2nd minor dim (" +
                    internalNodesParam.getMinorDimension2() + ") ! ");

        // minor dim 1 = nsite * 2
        int start = arrI * internalNodesParam.getMinorDimension1();
        for (int i = start; i < start + internalNodesParam.getMinorDimension1(); i++) {
            if (internalNodesParam.isDirty(i))
                return true;
        }
        return false;
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
    public List<String> getArguments() {
        return List.of();
//        return Collections.singletonList(internalNodesValues.getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
    public List<String> getConditions() {
        return List.of();
//        return internalNodesValues.getConditions();
    }

    @Override
    public void sample(State state, Random random) {}
}
