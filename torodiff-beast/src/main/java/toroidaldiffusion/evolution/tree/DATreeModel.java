package toroidaldiffusion.evolution.tree;

import beast.base.core.Description;
import beast.base.evolution.tree.TreeInterface;

@Description("Tree model for data augmentation tree likelihood")
public interface DATreeModel {

    TreeInterface getTree();

    // need Parameter.Base
    double[][] getTipValues();

    double[][] getInternalNodesValues();

    double[] getRootValues();

    int getSiteCount();

    int getLeafNodeCount();
}
