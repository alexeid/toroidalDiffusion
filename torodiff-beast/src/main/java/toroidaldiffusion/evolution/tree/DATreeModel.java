package toroidaldiffusion.evolution.tree;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;

@Description("Tree model for data augmentation tree likelihood")
public interface DATreeModel {

    TreeInterface getTree();

    RealParameter getTipsValuesParam();

    RealParameter getInternalNodesValuesParam();

    // need Parameter.Base
//    double[][] getTipsValues();

//    double[][] getInternalNodesValues();

//    double[] getRootValues();

    int getSiteCount();

    int getLeafNodeCount();

    double[] getNodeValue(Node node);
}
