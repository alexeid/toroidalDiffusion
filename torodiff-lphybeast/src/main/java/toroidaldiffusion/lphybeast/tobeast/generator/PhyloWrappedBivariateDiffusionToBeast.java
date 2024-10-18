package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.function.tree.InternalNodesID;
import lphy.core.model.Generator;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import toroidaldiffusion.DihedralAngleAlignment;
import toroidaldiffusion.Pair;
import toroidaldiffusion.PhyloWrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;
import toroidaldiffusion.lphybeast.tobeast.values.DihedralAnglesToBeast;
import toroidaldiffusion.operator.WrappedRandomWalkOperator;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static lphybeast.BEASTContext.getOperatorWeight;

public class PhyloWrappedBivariateDiffusionToBeast implements GeneratorToBEAST<PhyloWrappedBivariateDiffusion, toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion>{

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

        //sites
        Value<Double[][]> array = generator.getY();
        int site = (array.value().length)*2;

        //value sampled from generator, with internal nodes + tips
        DihedralAngleAlignment dihedralAngleAlignmentLphy = (DihedralAngleAlignment) dihedralAngleAlignmentValue.value();
        int taxa = dihedralAngleAlignmentLphy.pairs.length;

        /**
         * NO internalnodes is given from simulation
         * The internalnodes need to be sampled based on roots (?)
         */
        if (taxa == taxaNames.length()) {
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
        else if (taxa > taxaNames.length()) {

            tipValues = dihedralAngleAlignment;
            List<Double> internalNodesValues = getInternalNodes(dihedralAngleAlignmentValue);

            // Convert the list to a string with spaces between each value
            String internalNodesString = internalNodesValues.stream()
                    .map(String::valueOf)
                    .collect(Collectors.joining(" "));

            RealParameter internalNodes = new RealParameter(internalNodesString);

            //internalNodes.setInputValue("value", internalNodesString);
            internalNodes.setInputValue("minordimension", site);
            internalNodes.initAndValidate();

            context.addBEASTObject(internalNodes, dihedralAngleAlignmentValue);

            dihedralAngleTreeModel.setInputValue("tipValues", tipValues);
            dihedralAngleTreeModel.setInputValue("internalNodesValues", internalNodes);
        }

        else {
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
        removeNoIDTree(generator, context);

        dihedralAngleTreeModel.setInputValue("tree", timeTree);
        dihedralAngleTreeModel.initAndValidate();


        phyloWrappedBivariateDiffusion.setInputValue("daTreeModel", dihedralAngleTreeModel);

        phyloWrappedBivariateDiffusion.setInputValue("mu",
                context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.muParamName)));
        phyloWrappedBivariateDiffusion.setInputValue("sigma",
                context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.sigmaParamName)));
        phyloWrappedBivariateDiffusion.setInputValue("alpha",
                context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.alphaParamName)));

        WrappedRandomWalkOperator wrappedRandomWalkOperator = new WrappedRandomWalkOperator();
        RealParameter muParameter = context.getAsRealParameter(generator.getParams().get(PhyloWrappedBivariateDiffusion.muParamName));

        wrappedRandomWalkOperator.setInputValue("weight", getOperatorWeight(muParameter.getDimension() - 1));
        wrappedRandomWalkOperator.setInputValue("parameter", muParameter);
        wrappedRandomWalkOperator.setInputValue("windowSize", 0.1);
        //wrappedRandomWalkOperator.setInputValue("useGaussian", false);

        wrappedRandomWalkOperator.initAndValidate();
        wrappedRandomWalkOperator.setID("WrappedRandomWalkOperator." + muParameter.getID());

        context.addSkipOperator(muParameter);
        context.addExtraOperator(wrappedRandomWalkOperator);

        phyloWrappedBivariateDiffusion.initAndValidate();

        return phyloWrappedBivariateDiffusion;

    }

    public static List<Double> getInternalNodes(Value<DihedralAngleAlignment> dihedralAngleAlignmentValue) {

        List<Double> internalNodes = new ArrayList<>();

        DihedralAngleAlignment dihedralAngleAlignment = dihedralAngleAlignmentValue.value();

        String[] taxaNames = dihedralAngleAlignment.getTaxa().getTaxaNames();

        //pairsLength = tips + internalnodes
        int pairsLength = 2 * dihedralAngleAlignment.getTaxa().ntaxa() - 1;

        // Extract dihedral angles
        for (int i = taxaNames.length; i < pairsLength; i++) {
            for (int j = 0; j < dihedralAngleAlignment.nchar(); j++) {
                Pair pair = dihedralAngleAlignment.pairs[i][j];
                if (pair != null) {
                    // Add phi and psi values to the internalNodes list if they are non-null
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


    public void removeNoIDTree(PhyloWrappedBivariateDiffusion generator, BEASTContext context) {
        // Get the tree from the generator
        Value<TimeTree> tree = generator.getTree();

        // Get the BEAST object associated with the tree
        TreeInterface timeTree = (TreeInterface) context.getBEASTObject(tree);

        // if the tree's generator is not InternalNodesID (phi tree), then remove the tree
        if (!(tree.getGenerator() instanceof InternalNodesID)) {
            context.removeBEASTObject((BEASTInterface) timeTree);
        }
    }


    @Override
    public Class getGeneratorClass() {
        return PhyloWrappedBivariateDiffusion.class;
    }
}


