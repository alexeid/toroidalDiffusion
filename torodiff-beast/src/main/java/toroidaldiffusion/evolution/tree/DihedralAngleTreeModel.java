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

public class DihedralAngleTreeModel extends Distribution implements DATreeModel {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");

    // minor dimension 1 is sites, dim 2 is 2 (angles), dim = 2 * sites
    final public Input<RealParameter> tipValuesInput = new Input<>("tipValues", "");
    // assuming the root values are at the last
    final public Input<RealParameter> internalNodesValuesInput = new Input<>("internalNodesValues", "");

    // TODO use XOR rule to add y0 separately?
//    final public Input<RealParameter> rootValuesInput = new Input<>("rootValues", "y0, the value of [phi,psi] angle pairs for each carbon backbone bond of the molecule at the root of the phylogeny.");

    public final static int PAIR = 2;

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
    public double[][] getInternalNodesValues() {
        double[] values = internalNodesValuesInput.get().getDoubleValues();
        int nrows = this.getLeafNodeCount() - 1;
        int ncols = this.getSiteCount() * PAIR; // a pair
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
     *
     * @param node
     * @return  an array of pairs of values, dimension is 2 (angles) * N_sites
     */
    @Override
    public double[] getNodeValue(Node node) {
        if (node.isLeaf()) {
            String key = node.getID();

            Double[] values = getTipsValuesParam().getRowValues(key);
            return Stream.of(values).mapToDouble(Double::doubleValue).toArray();

        }
//            return getTipsValues()[node.getNr()];

        // TODO check: beast Nr starts ?
        return getInternalNodesValues()[node.getNr() - getLeafNodeCount()];
    }


    private double[][] convertTo2D(double[] flatArray, int nrows, int ncols) {
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

}
