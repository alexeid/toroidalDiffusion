package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import org.forester.applications.phylo2coloredgraphics;
import toroidaldiffusion.PhyloWrappedBivariateDiffusion;
import toroidaldiffusion.evolution.tree.DihedralAngleTreeModel;

import java.util.Arrays;

public class PhyloWrappedBivariateDiffusionToBeast implements GeneratorToBEAST<PhyloWrappedBivariateDiffusion, toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion>{


    @Override
    public toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion generatorToBEAST(PhyloWrappedBivariateDiffusion generator, BEASTInterface value, BEASTContext context) {

        toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion phyloWrappedBivariateDiffusion = new toroidaldiffusion.evolution.likelihood.PhyloWrappedBivariateDiffusion();

        //RealParameter dihedralAnglesAlignment = (RealParameter) value;
        //create DA model
        DihedralAngleTreeModel dihedralAngleTreeModel = new DihedralAngleTreeModel();

        //get the time tree
        Value<TimeTree> timeTreeValue = generator.getTree();
        TimeTree timeTree = timeTreeValue.value();
        dihedralAngleTreeModel.setInputValue("tree", timeTree);

        // Get the tip values (initial values for tips from the generator)
        Value<Double[][]> tipValues = generator.getY();
        Double[][] tipValuesArray = tipValues.value();

        RealParameter tipValuesParam = new RealParameter(Arrays.toString(tipValuesArray));
        dihedralAngleTreeModel.setInputValue("tipValues", tipValuesParam);

        //get internal nodes values
//        Double[][] internalNodesValuesArray = generator.getinternalNodes;
//        RealParameter internalNodesValuesParam = new RealParameter(internalNodesValuesArray);
//        dihedralAngleTreeModel.setInternalNodesValuesParam(internalNodesValuesParam);

        phyloWrappedBivariateDiffusion.setInputValue("daTreeModel", dihedralAngleTreeModel);


        return phyloWrappedBivariateDiffusion;
    }



    @Override
    public Class getGeneratorClass() {
        return PhyloWrappedBivariateDiffusion.class;
    }
}


