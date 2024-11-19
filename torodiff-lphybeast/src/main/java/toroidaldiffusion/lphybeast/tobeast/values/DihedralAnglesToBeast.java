package toroidaldiffusion.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
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

        /**
         * Get taxon names
         **/
        String[] taxaNames = dihedralAngleAlignment.getTaxa().getTaxaNames();
        StringBuilder builder = new StringBuilder();
        builder.append(taxaNames[0]);
        for (int i = 1; i < taxaNames.length; i++) {
            builder.append(" ");
            builder.append(taxaNames[i]);
        }


        int minordimension = dihedralAngleAlignment.nchar()*2;

        List<Double> angles = new ArrayList<>();

        /**
         * Extract dihedral angles tips
         */
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

        /**
         * DihedralAngles Only contains tips values
         */
        RealParameter dihedralAngles = new RealParameter();
        dihedralAngles.setInputValue("keys", builder.toString());
        dihedralAngles.setInputValue("value", angles);
        dihedralAngles.setInputValue("minordimension", minordimension);

        //TODO upper = "6.283"   2*pi
        dihedralAngles.setInputValue("upper", 6.283);

        dihedralAngles.initAndValidate();
        return dihedralAngles;
    }

    @Override
    public Class getValueClass() {
        return DihedralAngleAlignment.class;
    }
}

