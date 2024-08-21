package toroidaldiffusion.evolution.tree;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.State;

import java.util.List;
import java.util.Random;

public class DihedralAngleTreeModel extends Distribution {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");
    final public Input<Function> tipValuesInput = new Input<>("tipValues", "");
    final public Input<Function> internalNodesValuesInput = new Input<>("internalNodesValues", "");


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
