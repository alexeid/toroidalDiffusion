package toroidaldiffusion;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;

import java.util.Arrays;
import java.util.Objects;

/**
 *
 */
public class DihedralAngleAlignment implements TaxaCharacterMatrix<Pair> {

    // 1st[] is taxa, index is same order as Taxa
    // 2nd[] is the site 
    Pair[][] pairs;
    Taxa taxa;

    public DihedralAngleAlignment(Taxa taxa, int nchar) {
        this.taxa = taxa;
        this.pairs = new Pair[taxa.ntaxa()][nchar];
    }

    @Override
    public Pair getState(String taxonName, int column) {
        return pairs[taxa.indexOfTaxon(taxonName)][column];
    }

    @Override
    public void setState(int taxon, int position, Pair state) {
        pairs[taxon][position] = state;
    }

    @Override
    public Class getComponentType() {
        return Pair.class;
    }

    @Override
    public Taxa getTaxa() {
        return taxa;
    }

    @Override
    public Integer nchar() {
        return pairs[0].length;
    }

    @Override
    public String toJSON() {
        StringBuilder builder = new StringBuilder();
        builder.append("{\n");
        for (int i = 0; i < Objects.requireNonNull(taxa).ntaxa(); i++) {
            builder.append("  ").append(taxa.getTaxon(i));
            builder.append(" = ").append(Arrays.toString(pairs[i]));
//            if (i < n()-1)
            builder.append(",");
            builder.append("\n");
        }
        builder.append("  ntax = ").append(taxa.ntaxa());
        builder.append("\n").append("}");
        return builder.toString();
    }

    @Override
    public int getDimension() {
        return nchar()*taxa.getDimension();
    }

    public String toHTML() {
        return "";
    }
}
