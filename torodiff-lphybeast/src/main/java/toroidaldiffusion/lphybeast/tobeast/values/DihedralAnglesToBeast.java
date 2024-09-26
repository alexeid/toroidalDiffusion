package toroidaldiffusion.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.alignment.ContinuousCharacterData;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import toroidaldiffusion.DihedralAngleAlignment;
import toroidaldiffusion.Pair;

import java.util.ArrayList;
import java.util.List;


public class DihedralAnglesToBeast implements ValueToBEAST<DihedralAngleAlignment, RealParameter>{

    @Override
    public RealParameter valueToBEAST(Value<DihedralAngleAlignment> dihedralAngleAlignmentValue, BEASTContext context) {
        DihedralAngleAlignment dihedralAngleAlignment = dihedralAngleAlignmentValue.value();

        // Get taxa names
        String[] taxaNames = dihedralAngleAlignment.getTaxa().getTaxaNames();
        StringBuilder builder = new StringBuilder();
        builder.append(taxaNames[0]);
        for (int i = 1; i < taxaNames.length; i++) {
            builder.append(" ");
            builder.append(taxaNames[i]);
        }

        //Internal nodes
//        StringBuilder internalNodes = new StringBuilder();
//        int pairsLength = 2 * dihedralAngleAlignment.getTaxa().ntaxa() - 1;
//
//        for (int i = taxaNames.length; i < pairsLength; i++) {
//            for (int j = 0; j < dihedralAngleAlignment.nchar(); j++) {
//                Pair pair = dihedralAngleAlignment.getState(taxaNames[i], j);
//                if (pair != null) {
//                    if (pair.getPhi() != null) {
//                        internalNodes.append(pair.getPhi()).append("\t");
//                    }
//                    if (pair.getPsi() != null) {
//                        internalNodes.append(pair.getPsi()).append("\t");
//                    }
//                }
//            }
//            // Remove the trailing tab character and add a newline after each row
//            if (!internalNodes.isEmpty() && internalNodes.charAt(internalNodes.length() - 1) == '\t') {
//                internalNodes.deleteCharAt(internalNodes.length() - 1);
//            }
//            internalNodes.append("\n");
//
//            }

        List<Double> angles = new ArrayList<>();

        int minordimension = dihedralAngleAlignment.nchar()*2;

        // Extract dihedral angles
        for (String taxon : taxaNames) {
            for (int j = 0; j < dihedralAngleAlignment.nchar(); j++) {
                Pair pair = dihedralAngleAlignment.getState(taxon, j);
                if (pair != null) {
                    if (pair.getPhi() != null) {
                        angles.add(pair.getPhi());
                    }
                    if (pair.getPsi() != null) {
                        angles.add(pair.getPsi());
                    }
                }
            }
        }

        RealParameter dihedralAngles = new RealParameter();
        dihedralAngles.setInputValue("keys", builder.toString());
        dihedralAngles.setInputValue("value", angles);
        dihedralAngles.setInputValue("minordimension", minordimension);
        dihedralAngles.initAndValidate();

//        RealParameter internalNodesAngles = new RealParameter();
//        internalNodesAngles.setInputValue("values", Double.parseDouble(builder.toString()));
//        internalNodesAngles.setInputValue("minordimension", dihedralAngleAlignment.nchar());
//        internalNodesAngles.initAndValidate();

        return dihedralAngles;
    }

    @Override
    public Class getValueClass() {
        return DihedralAngleAlignment.class;
    }
}

