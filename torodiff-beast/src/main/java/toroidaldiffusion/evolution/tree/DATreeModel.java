package toroidaldiffusion.evolution.tree;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.evolution.tree.TreeInterface;

@Description("Tree model for data augmentation tree likelihood")
public interface DATreeModel {

    TreeInterface getTree();

    Function getTipValues();

    Function getInternalNodesValues();

}
