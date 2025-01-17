package toroidaldiffusion.evolution.tree;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;

import java.util.List;
import java.util.Random;
import java.util.stream.Stream;

/**
 * Use BEAST tree node Nr to map the internal node sequences (flatten in an 1d array parameter)
 * to the corresponding node during likelihood calculation.
 */
public class DihedralAngleTreeModel extends Distribution implements DATreeModel {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");

    // minor dimension 1 is sites, dim 2 is 2 (angles), dim = 2 * sites
    final public Input<RealParameter> tipValuesInput = new Input<>("tipValues", "");
    // assuming the root values are at the last
    final public Input<RealParameter> internalNodesValuesInput = new Input<>("internalNodesValues", "");

    // TODO use XOR rule to add y0 separately?
//    final public Input<RealParameter> rootValuesInput = new Input<>("rootValues", "y0, the value of [phi,psi] angle pairs for each carbon backbone bond of the molecule at the root of the phylogeny.");

    public final static int PAIR = 2;

    // index is node Nr, value is ID if node.getID() != null
    String[] nodeIDs;

    @Override
    public void initAndValidate() {
        // tree
        TreeInterface tree = treeInput.get();
        // data
        RealParameter tipValues = tipValuesInput.get();
        RealParameter internalNodesValues = internalNodesValuesInput.get();
        // validate dimensions
        final int leafNodeCount = tree.getLeafNodeCount();
        // getMinorDimension2 = Dimension() / this.minorDimension
        if (leafNodeCount != tipValues.getMinorDimension2())
            throw new IllegalArgumentException("The number of leaf nodes " + leafNodeCount +
                    " must be the same as the 2nd dimension of tip values " +
                    tipValues.getMinorDimension2() + " !");
        if (leafNodeCount - 1 != internalNodesValues.getMinorDimension2())
            throw new IllegalArgumentException("The number of leaf nodes - 1 " + (leafNodeCount-1) +
                    "must be the same as the 2nd dimension of internal nodes values ! " +
                    internalNodesValues.getMinorDimension2() + " !");

        List<String> keys = tipValues.getKeysList();
        if (keys.isEmpty() || keys.size() != leafNodeCount)
            throw new UnsupportedOperationException("Tip nodes value must use parameter having keys !");

        //TODO for getNodeValue(Node node)?
        nodeIDs = new String[tree.getNodeCount()];
        for (Node node : tree.getNodesAsArray()) {
            if (node.getID() != null)
                nodeIDs[node.getNr()] = node.getID();
            else
                nodeIDs[node.getNr()] = Integer.toString(node.getNr());
        }
    }

    public TreeInterface getTree() {
        return treeInput.get();
    }

    //TODO values array uses node.getNr() as index, need to check

    public RealParameter getTipsValuesParam() {
        return tipValuesInput.get();
    }

    public RealParameter getInternalNodesValuesParam() {
        return internalNodesValuesInput.get();
    }

    // TODO no keys?

    /**
     * @return  1st[] is the internal nodes index, assuming the last is the root.
     *          2nd[] are the pairs of angles flatten to 1d, ncols = nsite * 2;
     */
    public double[][] getInternalNodesValues() {
        // flatten 1d arr
        double[] values = internalNodesValuesInput.get().getDoubleValues();
        int nrows = this.getLeafNodeCount() - 1;
        int ncols = this.getSiteCount() * PAIR; // a pair
        // 1st[] is mapping to each node, 2nd[] is mapping to flatten 1d array for nsite*2 elements
        return convertTo2D(values, nrows, ncols);
    }

//    @Override
//    public double[] getRootValues() {
//        int len = getInternalNodesValues().length;
//        // assuming the root values are at the last
//        return getInternalNodesValues()[len-1];
//    }

    @Override
    public int getSiteCount() {
        return tipValuesInput.get().getMinorDimension1() / PAIR;
    }

    @Override
    public int getLeafNodeCount() {
        return getTree().getLeafNodeCount();
    }

    /**
     * Assuming the internal node sequences are arranged in the order of the BEAST node numbers Nr.
     * Additionally, the internal node sequences are stored in a 2D array,
     * with the number of rows corresponding to the total number of internal nodes.
     * @param node  a tip node or internal node.
     * @return  an array of pairs of values, its dimension is 2 (angles) * N_sites.
     */
    @Override
    public double[] getNodeValue(Node node) {
        if (node.isLeaf()) {
            // TODO key must start from 1, so node.getID() = node.getNr() + 1
            String key = node.getID();
            Double[] values = getTipsValuesParam().getRowValues(key);
            return Stream.of(values).mapToDouble(Double::doubleValue).toArray();
        }
        // beast Nr starts from 0, last is root
        int nr = node.getNr();
        int nTips = getLeafNodeCount();
        // TODO not use key ?
        return getInternalNodesValues()[nr - nTips]; // including root
    }

    private static double[][] convertTo2D(double[] flatArray, int nrows, int ncols) {
        if (flatArray.length != nrows * ncols)
            throw new IllegalArgumentException("Invalid length of tip flatArray ! " +
                    flatArray.length + " != " + nrows + " * " + ncols + " !");
        // Create the 2D array
        double[][] twoDArray = new double[nrows][ncols];
        // Fill the 2D array with elements from the 1D array
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                twoDArray[i][j] = flatArray[i * ncols + j];
            }
        }
        return twoDArray;
    }



    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {

    }

//    @Override
//    protected boolean requiresRecalculation() {
//        return treeInput.get().somethingIsDirty();
//    }

}
