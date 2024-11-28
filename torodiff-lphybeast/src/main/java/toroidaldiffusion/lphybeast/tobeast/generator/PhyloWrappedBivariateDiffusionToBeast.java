package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static toroidaldiffusion.PhyloWrappedBivariateDiffusion.y0RateParam;


public class PhyloWrappedBivariateDiffusionToBeast implements GeneratorToBEAST<PhyloWrappedBivariateDiffusion, toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion> {
    public static final double ANGLE_UPPER = 6.283; // 2 pi

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
        RealParameter tipValues;
        if (value instanceof RealParameter realParameter) {
            tipValues = realParameter;
        } else {
            throw new IllegalArgumentException("not a RealParameter");
        }
        dihedralAngleTreeModel.setInputValue("tipValues", tipValues);

        /**
         * use lphy DihedralAngleAlignment to determine if it contains the simulated internal node sequences,
         * if its number of sequences > ntaxa from tree, then it contains, otherwise nSeqs == ntaxa.
         */

        //extract Lphy (y) DA values which contains both tips + internalNodes
        Value dihedralAngleAlignmentValue = (Value) context.getGraphicalModelNode(tipValues);
        //Get the dihedralAngleAlignmentLphy values
        DihedralAngleAlignment dihedralAngleAlignmentLphy = (DihedralAngleAlignment) dihedralAngleAlignmentValue.value();

        //referenced taxa from tree (only tips)
        TimeTree tree = generator.getTree().value();
        // TODO do we need names ?
        Taxa taxaNames = tree.getTaxa();
        //Get internalnodes
        List<TimeTreeNode> internalNodes = tree.getInternalNodes();

        // minordimension = num of sites * 2
        int nsite = dihedralAngleAlignmentLphy.nchar();
        int minordimension =  nsite * 2;

        //value sampled from generator, with internal nodes + tips
        //nSeqs
        int nSeqs = dihedralAngleAlignmentLphy.pairs.length;

        RealParameter internalNodesSeqs;
        /**
         * NO internalnodes is given from simulation
         * The internalnodes need to be sampled based on roots (?)
         */
        if (nSeqs == taxaNames.length()) {
            Value y0V = generator.getParams().get(y0RateParam);

            // two options: 1. give root sequence, 2. give equilibrium frequency at the root
            if (y0V != null) {
                // root seqs
                Double[][] y0 = (Double[][]) y0V.value();
                // nsite rows, 2 cols
                Double[] flatArray = Arrays.stream(y0)
                        .flatMap(Arrays::stream)
                        .toArray(Double[]::new);

                int nINExclRoot = internalNodes.size() - 1;
                // exclude root
                String internalNodesString = "";
                if (nINExclRoot > 0) // there are other internal nodes
                    //TODO create internalNodes only contain root seq, rest uses 3.0 ?
                    internalNodesString = String.join(" ", Collections.nCopies(nINExclRoot * nsite, "3.14"));
                else if (nINExclRoot < 0)
                    throw new IllegalStateException("Incorrect number of internal nodes : " + internalNodes.size());
                // root seq always at the last
                String rootSeqStr = String.join(" ", Arrays.stream(flatArray)
                        .map(String::valueOf) // Convert Double to String
                        .toArray(String[]::new) );
                internalNodesString += rootSeqStr;

                //TODO err

                internalNodesSeqs = new RealParameter(internalNodesString);
                //internalNodes.setInputValue("value", internalNodesString);
                internalNodesSeqs.setInputValue("minordimension", minordimension );
                internalNodesSeqs.initAndValidate();

                //TODO add WrpRanOP here to sample internal node seq


            } else
                throw new UnsupportedOperationException("Giving equilibrium frequency at the root is not supported yet !");

        } else if (nSeqs > taxaNames.length()) {
            //Get nodesID
            String[] internalNodesID = getID(internalNodes);

            // Internal nodes sequences are given from simulation
            internalNodesSeqs = getInternalNodesParam(dihedralAngleAlignmentValue, internalNodesID, internalNodes, minordimension);
            context.addBEASTObject(internalNodesSeqs, dihedralAngleAlignmentValue);
        } else {
            throw new IllegalStateException("Unexpected condition: 'taxa' is less than expected or invalid.");
        }

        dihedralAngleTreeModel.setInputValue("internalNodesValues", internalNodesSeqs);

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

        //*** https://github.com/LinguaPhylo/LPhyBeast/issues/172 **//
        if (timeTree instanceof TreeParser treeParser)
            // IsLabelledNewick must = false, otherwise internal node index will be messed up
            treeParser.setInputValue("IsLabelledNewick", false);

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
        muParameter.setInputValue("upper", ANGLE_UPPER);
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

    public static RealParameter getInternalNodesParam(Value dihedralAngleAlignmentValue, String[] internalNodesID, List<TimeTreeNode> internalNodesList, int minordimension) {
        List<Double> internalNodesValues = getInternalNodes(dihedralAngleAlignmentValue, internalNodesID, internalNodesList);
        // Convert list to a string with spaces between each value for input values
        String internalNodesString = internalNodesValues.stream()
                .map(String::valueOf)
                .collect(Collectors.joining(" "));
        // Convert string list to string for keys
        String internalNodeskeys = String.join(" ", internalNodesID);

        return createInternalNodeValuesParameter(minordimension, internalNodesString, internalNodeskeys);
    }

    private static RealParameter createInternalNodeValuesParameter(int minordimension, String internalNodesString,
                                                                   String internalNodeskeys) {
        RealParameter internalNodes = new RealParameter(internalNodesString);
        //internalNodes.setInputValue("value", internalNodesString);
        internalNodes.setInputValue("minordimension", minordimension); // a pair
        if (internalNodeskeys != null)
            internalNodes.setInputValue("keys", internalNodeskeys);
        //TODO upper = "6.283"   2*pi
        internalNodes.setInputValue("upper", ANGLE_UPPER);
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




