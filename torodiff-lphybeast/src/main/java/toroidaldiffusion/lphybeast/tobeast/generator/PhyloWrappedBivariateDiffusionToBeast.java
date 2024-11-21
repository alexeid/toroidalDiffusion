package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.operators.BICESPSTreeOperatorStractegy;
import lphybeast.tobeast.operators.TreeOperatorStrategy;
import toroidaldiffusion.DihedralAngleAlignment;
import toroidaldiffusion.Pair;
import toroidaldiffusion.PhyloWrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;



public class PhyloWrappedBivariateDiffusionToBeast implements GeneratorToBEAST<PhyloWrappedBivariateDiffusion, toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion> {

    //without internalnodes sequence - cannot calculate likelihood
    //the internalnodes involve the root initialised angles

    @Override
    public toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion generatorToBEAST(PhyloWrappedBivariateDiffusion generator, BEASTInterface value, BEASTContext context) {

        toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion phyloWrappedBivariateDiffusion = new toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion();

        /**
         * Create DA model
         */
        DihedralAngleTreeModel dihedralAngleTreeModel = new DihedralAngleTreeModel();

        /**
         * dihedralAngleAlignment only contains tips
         */
        RealParameter dihedralAngleAlignment;
        if (value instanceof RealParameter realParameter) {
            dihedralAngleAlignment = realParameter;
        } else {
            throw new IllegalArgumentException("not a RealParameter");
        }

        RealParameter tipValues;

        //extract Lphy (y) DA values which contains both tips + internalNodes
        Value dihedralAngleAlignmentValue = (Value) context.getGraphicalModelNode(dihedralAngleAlignment);

        //referenced taxa from tree (only tips)
        TimeTree tree = generator.getTree().value();
        Taxa taxaNames = tree.getTaxa();

        //Get internalnodes
        List<TimeTreeNode> internalNodes = tree.getInternalNodes();

        //Get the dihedralAngleAlignmentLphy values
        DihedralAngleAlignment dihedralAngleAlignmentLphy = (DihedralAngleAlignment) dihedralAngleAlignmentValue.value();

        //Get nodesID
        String[] internalNodesID = getID(internalNodes);

        //sites
        Value<Double[][]> array = generator.getY();
        int site = (array.value().length) * 2;

        //value sampled from generator, with internal nodes + tips
        //nSeqs
        int nSeqs = dihedralAngleAlignmentLphy.pairs.length;

        /**
         * NO internalnodes is given from simulation
         * The internalnodes need to be sampled based on roots (?)
         */
        if (nSeqs == taxaNames.length()) {
//            tipValues = dihedralAngleAlignment;
//            System.out.println("no internalnodes given");
//            //todo: if no internalnodes sequences given, internalnodes sequences need to be sampled based on root sequence
//            Value<Double[][]> rootSequence = generator.getY();
//            internalNodes = rootSequence.sampling...

            throw new IllegalStateException("Unexpected condition: 'taxa' is less than expected or invalid.");

        }

        /**
         * Internalnodes is given from simulation
         */
        else if (nSeqs > taxaNames.length()) {

            tipValues = dihedralAngleAlignment;

            RealParameter internalNodesValues = setInternalNodesRP(dihedralAngleAlignmentValue, internalNodesID, internalNodes, site);
            //TODO upper = "6.283"   2*pi
            internalNodesValues.setInputValue("upper", 6.283);

            context.addBEASTObject(internalNodesValues, dihedralAngleAlignmentValue);

            dihedralAngleTreeModel.setInputValue("tipValues", tipValues);
            dihedralAngleTreeModel.setInputValue("internalNodesValues", internalNodesValues);
        } else {
            throw new IllegalStateException("Unexpected condition: 'taxa' is less than expected or invalid.");
        }

        //set inputs for DAtree model
//        if (internalNodes == null) {
//            dihedralAngleTreeModel.setInputValue("tipValues", dihedralAngleAlignment);
//            dihedralAngleTreeModel.setInputValue("internalNodes", internalNodes);
//            //todo: add sampled internalnodes
//            // dihedralAngleTreeModel.setInputValue("internalNodes", internalNodes);
//        }
//        else {
//            dihedralAngleTreeModel.setInputValue("tipValues", tipValues);
//            dihedralAngleTreeModel.setInputValue("internalNodes", internalNodes);
//        }

        /**
         * Get the Timetree
         */
        TreeInterface timeTree = (TreeInterface) context.getBEASTObject(generator.getTree());
        //removeNoIDTree(generator, context);

        dihedralAngleTreeModel.setInputValue("tree", timeTree);
        dihedralAngleTreeModel.initAndValidate();

        //*** TODO for dev ***
        Tree stateNodeTree = (Tree) timeTree;
        context.addSkipOperator(stateNodeTree);
        // fix tree topology
        Operator rootAgeScale = TreeOperatorStrategy.createRootHeightOperator(stateNodeTree, context);
        Operator treeFlex = BICESPSTreeOperatorStractegy.createBICEPSTreeFlex(stateNodeTree, context);
        context.addExtraOperator(rootAgeScale);
        context.addExtraOperator(treeFlex);
        //*** end ***

        phyloWrappedBivariateDiffusion.setInputValue("daTreeModel", dihedralAngleTreeModel);

        RealParameter muParameter = context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.muParamName));
        //TODO upper = "6.283"   2*pi
        muParameter.setInputValue("upper", 6.283);
        phyloWrappedBivariateDiffusion.setInputValue("mu", muParameter);
        phyloWrappedBivariateDiffusion.setInputValue("sigma",
                context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.sigmaParamName)));


        phyloWrappedBivariateDiffusion.setInputValue("drift", context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.DRIFT_PARAM)));
        phyloWrappedBivariateDiffusion.setInputValue("driftCorr", context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.DRIFT_CORR_PARAM)));

// TODO rm, this op is used for sampling internal node sequences
//        WrappedRandomWalkOperator wrappedRandomWalkOperator = new WrappedRandomWalkOperator();
//        RealParameter muParameter = context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.muParamName));
//
//        wrappedRandomWalkOperator.setInputValue("weight", getOperatorWeight(muParameter.getDimension()));
//        wrappedRandomWalkOperator.setInputValue("parameter", muParameter);
//        wrappedRandomWalkOperator.setInputValue("windowSize", 0.1);
//
//        wrappedRandomWalkOperator.initAndValidate();
//        wrappedRandomWalkOperator.setID("WrappedRandomWalkOperator." + muParameter.getID());

//        context.addSkipOperator(muParameter);
//        context.addExtraOperator(wrappedRandomWalkOperator);

        phyloWrappedBivariateDiffusion.initAndValidate();

        context.addExtraLoggable(phyloWrappedBivariateDiffusion);

        return phyloWrappedBivariateDiffusion;

    }

    public static RealParameter setInternalNodesRP(Value dihedralAngleAlignmentValue, String[] internalNodesID, List<TimeTreeNode> internalNodesList, int site) {
        List<Double> internalNodesValues = getInternalNodes(dihedralAngleAlignmentValue, internalNodesID, internalNodesList);
        // Convert list to a string with spaces between each value for input values
        String internalNodesString = internalNodesValues.stream()
                .map(String::valueOf)
                .collect(Collectors.joining(" "));
        // Convert string list to string for keys
        String internalNodeskeys = String.join(" ", internalNodesID);

        RealParameter internalNodes = new RealParameter(internalNodesString);
        //internalNodes.setInputValue("value", internalNodesString);
        internalNodes.setInputValue("minordimension", site);
        internalNodes.setInputValue("keys", internalNodeskeys);
        internalNodes.initAndValidate();
        return internalNodes;
    }
    

    public String[] getID(List<TimeTreeNode> nodes) {
        // Return empty array if the input list is null or empty
        if (nodes == null || nodes.isEmpty()) {
            return new String[0];
        }

        // Create an array of the same size as the list
        String[] ids = new String[nodes.size()];

        for (int i = 0; i < nodes.size(); i++) {
            TimeTreeNode node = nodes.get(i);
            String id = node.getId();
            ids[i] = (id != null) ? id : "";
        }

        return ids;
    }
    

    public static List<Double> getInternalNodes(Value<DihedralAngleAlignment> dihedralAngleAlignmentValue, String[] nodeID, List<TimeTreeNode> nodes) {
        DihedralAngleAlignment dihedralAngleAlignment = dihedralAngleAlignmentValue.value();

        Taxa taxa = dihedralAngleAlignment.getTaxa();

        List<Double> internalNodes = new ArrayList<>();

        for (int i = 0; i < nodes.size(); i++) {
            int nodeIndex;
            nodeIndex = dihedralAngleAlignment.indexOfSequence(nodeID[i], taxa, nodes);
            for (int j = 0; j < dihedralAngleAlignment.nchar(); j++) {
                Pair pair = dihedralAngleAlignment.getState(nodeIndex, j);
                if (pair != null) {
                    if (pair.getPhi() != null) {
                        internalNodes.add(pair.getPhi());
                    }
                    if (pair.getPsi() != null) {
                        internalNodes.add(pair.getPsi());
                    }
                }

            }
        }

        return internalNodes;
    }



    @Override
    public Class getGeneratorClass() {
        return PhyloWrappedBivariateDiffusion.class;
    }
}




