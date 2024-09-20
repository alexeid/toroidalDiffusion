package tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import toroidaldiffusion.DihedralAngleAlignment;
import toroidaldiffusion.Pair;


public class DihedralAnglesToBeast implements ValueToBEAST<DihedralAngleAlignment, RealParameter>{

    @Override
    public RealParameter valueToBEAST(Value<DihedralAngleAlignment> dihedralAngleAlignmentValue, BEASTContext context) {
        DihedralAngleAlignment dihedralAngleAlignment = dihedralAngleAlignmentValue.value();

        // Get taxa names
        String[] taxaNames = dihedralAngleAlignment.getTaxa().getTaxaNames();

        // DihedralAngles alignment
        StringBuilder builder = new StringBuilder();
        for (String taxon : taxaNames) {
            for (int j = 0; j < dihedralAngleAlignment.nchar(); j++) {
                // Get the dihedral angles pair (phi, psi) for the current taxon and position
                Pair pair = dihedralAngleAlignment.getState(taxon, j);
                if (pair != null) {
                    if (pair.getPhi() != null) {
                        builder.append(pair.getPhi()).append("\t");
                    }
                    if (pair.getPsi() != null) {
                        builder.append(pair.getPsi()).append("\t");
                    }
                }
            }
            // Remove the trailing tab character and add a newline after each row
            if (!builder.isEmpty() && builder.charAt(builder.length() - 1) == '\t') {
                builder.deleteCharAt(builder.length() - 1);
            }
            builder.append("\n");
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

        RealParameter dihedralAngles = new RealParameter();
        dihedralAngles.setInputValue("keys", taxaNames);
        dihedralAngles.setInputValue("values", Double.parseDouble(builder.toString()));
        dihedralAngles.setInputValue("minordimension", dihedralAngleAlignment.nchar());
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
