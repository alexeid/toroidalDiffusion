package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import lphy.core.parser.LPhyListenerImpl;
import lphy.core.parser.LPhyParserDictionary;
import lphy.core.parser.REPL;
import org.junit.jupiter.api.Test;
import toroidaldiffusion.DihedralAngleAlignment;
import toroidaldiffusion.Pair;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;

class PhyloWrappedBivariateDiffusionToBeastTest {

    @Test
    public void getRealParameter() {
        LPhyListenerImpl parser = new LPhyListenerImpl(new REPL());
        String Lphymodel = "model {\n  λ ~ LogNormal(meanlog=3.0, sdlog=0.1);\n  " +
                "ψ ~ Yule(lambda=λ, n=10);\n  " +
                "ψ2 = setInternalNodesID(ψ);\n  " +
                "mu = [2.0, 2.5];\n  " +
                "alpha = [1.0, 1.0, 0.5];\n  " +
                "σ = [0.5, 0.75];\n  " +
                "y0 = [[3.0, 3.0], [3.0, 3.0]];\n  " +
                "y ~ PhyloWrappedBivariateDiffusion(mu=mu, sigma=σ, alpha=alpha, y0=y0, tree=ψ2);\n}";

        parser.parse(Lphymodel);

        LPhyParserDictionary parserDictionary = parser.getParserDictionary();

        Map<String, Value<?>> modelDict = parserDictionary.getModelDictionary();
        Set<Value> modelValueSet = parserDictionary.getModelValues();

        Value dihedralAngleValue = modelDict.get("y");
        DihedralAngleAlignment dihedralAngleAlignment = (DihedralAngleAlignment) dihedralAngleValue.value();
        Value treeValue = modelDict.get("ψ2");
        TimeTree tree = (TimeTree) treeValue.value();

        //get internalnodes
        List<TimeTreeNode> internalNodes = tree.getInternalNodes();

        //get id stringlist
        String[] ids = new String[internalNodes.size()];

        for (int i = 0; i < internalNodes.size(); i++) {
            TimeTreeNode node = internalNodes.get(i);
            String id = node.getId();
            ids[i] = (id != null) ? id : "";
        }

        int site = dihedralAngleAlignment.nchar() * 2;

        RealParameter internalNodesValues = PhyloWrappedBivariateDiffusionToBeast.getInternalNodesParam(dihedralAngleValue, ids, internalNodes, site);

        //check for the ids
        String expectedKeys = String.join(" ", ids);

        //get expected values
        String[] taxaNames = dihedralAngleAlignment.getTaxa().getTaxaNames();
        List<Double> expectedValues = new ArrayList<>();

        int pairsLength = 2 * dihedralAngleAlignment.getTaxa().ntaxa() - 1;
        for (int i = taxaNames.length; i < pairsLength; i++) {
            for (int j = 0; j < dihedralAngleAlignment.nchar(); j++) {
                Pair pair = dihedralAngleAlignment.getState(i, j);
                if (pair != null) {
                    System.out.println("Pair found at (" + i + ", " + j + ")");
                    if (pair.getPhi() != null) {
                        expectedValues.add(pair.getPhi());
                        System.out.println("Phi value: " + pair.getPhi());
                    }
                    if (pair.getPsi() != null) {
                        expectedValues.add(pair.getPsi());
                        System.out.println("Psi value: " + pair.getPsi());
                    }
                }
            }
        }

                for (int i2 = 0; i2 < expectedValues.size(); i2++) {
                    assertEquals(expectedValues.get(i2), internalNodesValues.getValue(i2), "Value mismatch at index " + i2);
                }

                assertEquals(site, internalNodesValues.getInputValue("minordimension"), "Incorrect site dimension");
                assertEquals(expectedKeys, internalNodesValues.getInputValue("keys"), "Incorrect keys for internal nodes");


            }

        }
