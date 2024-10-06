package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import toroidaldiffusion.PhyloWrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;

public class PhyloWrappedBivariateDiffusionToBeast implements GeneratorToBEAST<PhyloWrappedBivariateDiffusion, toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion>{

    //without internalnodes sequence - cannot calculate likelihood
    //the internalnodes involve the root initialised angles

    @Override
    public toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion generatorToBEAST(PhyloWrappedBivariateDiffusion generator, BEASTInterface value, BEASTContext context) {

        toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion phyloWrappedBivariateDiffusion = new toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion();

        //create DA model
        DihedralAngleTreeModel dihedralAngleTreeModel = new DihedralAngleTreeModel();


        RealParameter dihedralAngleAlignment;
        if (value instanceof RealParameter realParameter) {
            dihedralAngleAlignment = realParameter;
        } else {
            throw new IllegalArgumentException("not a RealParameter");
        }

        RealParameter tipValues;
        RealParameter internalNodes = null;

        Value dihedralAngleAlignmentValue = (Value) context.getGraphicalModelNode(dihedralAngleAlignment);

        //reference taxa length
        TimeTree tree = generator.getTree().value();
        Taxa taxaNames = tree.getTaxa();
        //String taxaNames = generator.getTree().getId();

        Value<Double[][]> array = generator.getY();
        int site = array.value().length;
        int taxa = dihedralAngleAlignment.getDimension()/site/2;

        if (dihedralAngleAlignment == dihedralAngleAlignmentValue.value()) {
            tipValues = dihedralAngleAlignment;

            System.out.println("no internalnodes given");
            //todo: if no internalnodes sequences given, internalnodes sequences need to be sampled based on root sequence
//            Value<Double[][]> rootSequence = generator.getY();
//            internalNodes = rootSequence.sampling...
        }
        else if (taxa > taxaNames.length()) {
            tipValues = new RealParameter();
            for (int i = 0; i < (taxaNames.length() * site * 2); i++) {
                double angle = dihedralAngleAlignment.getValue(i);
                tipValues.setValue(i, angle);
            }

            internalNodes = new RealParameter();
            for (int i = taxaNames.length() * site * 2; i < dihedralAngleAlignment.getDimension(); i++) {
                double angle = dihedralAngleAlignment.getValue(i);
                internalNodes.setValue(i, angle);
            }

            context.removeBEASTObject((BEASTInterface) dihedralAngleAlignmentValue);
            context.addBEASTObject(tipValues, dihedralAngleAlignmentValue);
            context.addBEASTObject(internalNodes, dihedralAngleAlignmentValue);
        }
        else {
            throw new IllegalStateException("Unexpected condition: 'taxa' is less than expected or invalid.");
        }

        //set inputs for DAtree model
        if (internalNodes == null) {
            dihedralAngleTreeModel.setInputValue("tipValues", dihedralAngleAlignment);
        //todo: add sampled internalnodes
            // dihedralAngleTreeModel.setInputValue("internalNodes", internalNodes);
        }
        else {
            dihedralAngleTreeModel.setInputValue("tipValues", tipValues);
            dihedralAngleTreeModel.setInputValue("internalNodes", internalNodes);
        }

        //get the time tree
        TreeInterface timeTree = (TreeInterface) context.getBEASTObject(generator.getTree());
        dihedralAngleTreeModel.setInputValue("tree", timeTree);
        dihedralAngleTreeModel.initAndValidate();


        phyloWrappedBivariateDiffusion.setInputValue("daTreeModel", dihedralAngleTreeModel);
        phyloWrappedBivariateDiffusion.setInputValue("mu", context.getAsRealParameter(generator.getParams().get("mu")));
        phyloWrappedBivariateDiffusion.setInputValue("sigma", context.getAsRealParameter(generator.getParams().get("sigma")));
        phyloWrappedBivariateDiffusion.setInputValue("alpha", context.getAsRealParameter(generator.getParams().get("alpha")));

        phyloWrappedBivariateDiffusion.initAndValidate();

        return phyloWrappedBivariateDiffusion;
    }



    @Override
    public Class getGeneratorClass() {
        return PhyloWrappedBivariateDiffusion.class;
    }
}


