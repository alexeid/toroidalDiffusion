package toroidaldiffusion.evolution.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
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
     * internal nodes values
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
