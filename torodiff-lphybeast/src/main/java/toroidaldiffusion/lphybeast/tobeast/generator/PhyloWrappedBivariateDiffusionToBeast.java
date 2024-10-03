package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.GraphicalModelNode;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import org.forester.applications.phylo2coloredgraphics;
import toroidaldiffusion.DihedralAngleAlignment;
import toroidaldiffusion.PhyloWrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;
import toroidaldiffusion.lphybeast.tobeast.values.DihedralAnglesToBeast;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PhyloWrappedBivariateDiffusionToBeast implements GeneratorToBEAST<PhyloWrappedBivariateDiffusion, toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion>{


    @Override
    public toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion generatorToBEAST(PhyloWrappedBivariateDiffusion generator, BEASTInterface value, BEASTContext context) {

        final RealParameter dihedralAngleAlignment;
        if (value instanceof RealParameter realParameter) {
            dihedralAngleAlignment = realParameter;
        }
        else {
            throw new IllegalArgumentException ("not a RealParameter");
        }


        toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion phyloWrappedBivariateDiffusion = new toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion();

        //create DA model
        DihedralAngleTreeModel dihedralAngleTreeModel = new DihedralAngleTreeModel();

        //get the time tree
        //Value<TimeTree> timeTreeValue = generator.getTree();
        TreeInterface timeTree = (TreeInterface) context.getBEASTObject(generator.getTree());
        dihedralAngleTreeModel.setInputValue("tree", timeTree);

        // Get the tip values (initial values for tips from the generator)
       //y0 = root values, not tips value, need to figure out tipsvalues
        //RealParameter tipValuesParam = context.getAsRealParameter(DihedralAnglesToBeast.);

        //dihedralAngleTreeModel.setInputValue("tipValues", dihedralAngleAlignmentValue);

        //get internal nodes values
//        Double[][] internalNodesValuesArray = generator.getinternalNodes;
//        RealParameter internalNodesValuesParam = new RealParameter(internalNodesValuesArray);
//        dihedralAngleTreeModel.setInternalNodesValuesParam(internalNodesValuesParam);

        // int pairsLength = 2 * dihedralAngleAlignment.getTaxa().ntaxa() - 1;
        // get the num of taxa

        //get the num of taxa from values
        //num of sites = length of y0 array
        Value<Double[][]> array = generator.getY();
        int site = array.value().length;
        int taxa = dihedralAngleAlignment.getDimension()/site/2;

        //get the actual num of taxa
        DihedralAngleAlignment dihedralAngleAlignmentValue = (DihedralAngleAlignment) value; //TODO: can we cast realparameter -> DihedralAngleAlignment?
        String [] taxaNames = dihedralAngleAlignmentValue.getTaxa().getTaxaNames();

        //if the number of taxa equals to the taxaNames length, tipsValues == dihedralAngleAlignment stored values
        if (taxa == taxaNames.length) {
            dihedralAngleTreeModel.setInputValue("tipValues", dihedralAngleAlignment);
        }
//        else if (taxa > taxaNames.length) {
//            List<Double> tipsValues = new ArrayList<>();
//            //num of tips angles = taxon * site * 2 (pairs)
//            for (int i = 0; i < taxaNames.length*site*2 - 1; i++) {
//                double angle = dihedralAngleAlignment.getValue(i);
//                tipsValues.add(angle);
//            }
//            List<Double> internalNodes = new ArrayList<>();
//            //values beyond the tip values are internal node values
//            for (int i = taxaNames.length*site*2; i < (dihedralAngleAlignment.getDimension() - 1); i++) {
//                double angle = dihedralAngleAlignment.getValue(i);
//                internalNodes.add(angle);
//            }
//
//            dihedralAngleTreeModel.setInputValue("tipValues", context.getAsRealParameter((Value) tipsValues));
//            dihedralAngleTreeModel.setInputValue("internalNodesAngles", context.getAsRealParameter((Value) internalNodes));
//
//            //todo: how to keep dihedralAngleAlignment as input but change its value to keep only tipsangles
//            //todo: context.addBEASTObject() required inputs are the arguments of a function or distribution or the function/distribution
//        }

        //if the number of taxa is greater than taxaNames length, then:
        //1. from index 0 to index to (taxaNames.length * site * 2) == tips values
        //2. from index (taxaNames.length * site * 2) to index realparameter dimension == internal nodes values
        else if (taxa > taxaNames.length) {
            for (int i = 0; i < (taxaNames.length * site * 2); i++) {
                double angle = dihedralAngleAlignment.getValue(i);
                dihedralAngleAlignment.setValue(i, angle);
            }

            RealParameter internalnodes = new RealParameter();
            for (int i = taxaNames.length * site * 2; i < dihedralAngleAlignment.getDimension(); i++) {
                double angle = dihedralAngleAlignment.getValue(i);
                internalnodes.setValue(i, angle);
            }
            dihedralAngleTreeModel.setInputValue("tipValues", dihedralAngleAlignment);
            context.addBEASTObject(dihedralAngleTreeModel, (GraphicalModelNode) internalnodes);
            //todo: context.addBEASTObject() required inputs are the arguments of a function or distribution or the function/distribution
        }

        //if is neither of the conditions above, the dihedralAngleAlignment fails
        else {
            throw new IllegalStateException("Unexpected condition: 'taxa' is less than expected or invalid.");
        }

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


