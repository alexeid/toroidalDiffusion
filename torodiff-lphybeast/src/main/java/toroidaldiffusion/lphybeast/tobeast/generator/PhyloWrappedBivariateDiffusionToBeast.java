package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import lphy.base.distribution.Uniform;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Generator;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import lphybeast.tobeast.operators.BICESPSTreeOperatorStractegy;
import lphybeast.tobeast.operators.TreeOperatorStrategy;
import toroidaldiffusion.DihedralAngleAlignment;
import toroidaldiffusion.Pair;
import toroidaldiffusion.PhyloWrappedBivariateDiffusion;
import toroidaldiffusion.WrappedNormalConst;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;
import toroidaldiffusion.operator.WrappedRandomWalkOperator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static lphybeast.BEASTContext.getOperatorWeight;
import static toroidaldiffusion.WrappedNormalConst.MAX_ANGLE_VALUE;


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
//TODO root (all internal nodes?) sequences are sampled from equilibrium distribution
                throw new UnsupportedOperationException("Not supported yet !");

        } else if (nSeqs > taxaNames.length()) {
            //Get nodesID
            String[] internalNodesID = getID(internalNodes);

            boolean sampleInterNodeSeq = true; //TODO how get this flag from lphy?
            if (sampleInterNodeSeq) {
                // TODO cannot handle keys if sampling internal node seqs
                internalNodesSeqs = getInternalNodesParam(dihedralAngleAlignmentValue, internalNodesID,
                        internalNodes, minordimension, false);

                context.addStateNode(internalNodesSeqs, dihedralAngleAlignmentValue, false);
                //TODO add WrappedRandomWalkOperator to sample internal node seq
                WrappedRandomWalkOperator wrappedRWOp = new WrappedRandomWalkOperator();
                wrappedRWOp.initByName("parameter", internalNodesSeqs, "windowSize", 0.1,
                        "weight", getOperatorWeight(internalNodesSeqs.getDimension() - 1));
                // TODO tree id?
                wrappedRWOp.setID("" + "InternalNodeSeqs.WrappedRandomWalk");
                context.addExtraOperator(wrappedRWOp);

                //TODO rm logging internal node seqs
                context.addSkipLoggable(internalNodesSeqs);

            } else
                // Internal nodes sequences are given from simulation
                internalNodesSeqs = getInternalNodesParam(dihedralAngleAlignmentValue, internalNodesID,
                        internalNodes, minordimension, true);

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

        Map<String, Value> paramMap = generator.getParams();
        RealParameter muParameter = context.getAsRealParameter(paramMap.get(WrappedNormalConst.muParamName));
        //TODO upper = "6.283"   2*pi
        muParameter.setInputValue("upper", MAX_ANGLE_VALUE);
        phyloWrappedBivariateDiffusion.setInputValue("mu", muParameter);
        phyloWrappedBivariateDiffusion.setInputValue("sigma",
                context.getAsRealParameter(paramMap.get(WrappedNormalConst.sigmaParamName)));

        phyloWrappedBivariateDiffusion.setInputValue("drift", context.getAsRealParameter(paramMap.get(WrappedNormalConst.DRIFT_PARAM)));

        Value driftCorrValue = paramMap.get(WrappedNormalConst.DRIFT_CORR_PARAM);
        RealParameter driftCorrParam = context.getAsRealParameter(driftCorrValue);
        Generator genDC = driftCorrValue.getGenerator();
        //TODO why beast Uniform prior not restrict the bound ?
        if (genDC instanceof Uniform uniformLPhy) {
            driftCorrParam.initByName("lower", ValueUtils.doubleValue(uniformLPhy.getLower()),
                    "upper", ValueUtils.doubleValue(uniformLPhy.getUpper()));
        } else
            driftCorrParam.initByName("lower", -1.0, "upper", 1.0);
        phyloWrappedBivariateDiffusion.setInputValue("driftCorr", driftCorrParam);

// TODO rm, this op is used for sampling internal node sequences
//        WrappedRandomWalkOperator wrappedRandomWalkOperator = new WrappedRandomWalkOperator();
//        RealParameter muParameter = context.getAsRealParameter(paramMap.get(PhyloWrappedBivariateDiffusion.muParamName));
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

    public static RealParameter getInternalNodesParam(
            Value dihedralAngleAlignmentValue, String[] internalNodesID,
            List<TimeTreeNode> internalNodesList, int minordimension, boolean useKey) {

        List<Double> internalNodesValues = getInternalNodes(dihedralAngleAlignmentValue, internalNodesID, internalNodesList);
        // Convert list to a string with spaces between each value for input values
        String internalNodesString = internalNodesValues.stream()
                .map(String::valueOf)
                .collect(Collectors.joining(" "));

        RealParameter internalNodes = new RealParameter(internalNodesString);
        //internalNodes.setInputValue("value", internalNodesString);
        // minordimension = num of sites * 2
        internalNodes.setInputValue("minordimension", minordimension);
        if (useKey) {
            // Convert string list to string for keys
            String internalNodeskeys = String.join(" ", internalNodesID);
            internalNodes.setInputValue("keys", internalNodeskeys);
        }

        //TODO upper = "6.283"   2*pi
        internalNodes.setInputValue("upper", MAX_ANGLE_VALUE);
        internalNodes.initAndValidate();
        // TODO unique?
        internalNodes.setID("internalNodes");
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




