package toroidaldiffusion;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.core.logger.TextFileFormatted;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 *
 */
public class DihedralAngleAlignment implements TaxaCharacterMatrix<Pair>, TextFileFormatted {

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
            builder.append("  ").append(taxa.getTaxon(i)); //getname not tostring
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
        StringBuilder builder = new StringBuilder();

        builder.append("<html>\n");
        builder.append("<body>\n");
        builder.append("<table border='1'>\n");

        // Add table headers
        builder.append("<tr>");
        builder.append("<th>Taxon</th>");
        for (int j = 0; j < nchar(); j++) {
            builder.append("<th>Site ").append(j + 1).append("</th>");
        }
        builder.append("</tr>\n");

        // Add table rows for each taxon
//            for (int i = 0; i < taxa.ntaxa(); i++) {
//                builder.append("<tr>");
//                builder.append("<td>").append(taxa.getTaxon(i)).append("</td>");
//                for (int j = 0; j < nchar(); j++) {
//                    builder.append("<td>").append(pairs[i][j] != null ? pairs[i][j].toString() : "").append("</td>");
//                }
//                builder.append("</tr>\n");
//            }

        String[] taxaNames = taxa.getTaxaNames();
        for (int i = 0; i < taxaNames.length; i++) {
            builder.append("<tr>");
            builder.append("<td>").append(taxaNames[i]).append("</td>");
            for (int j = 0; j < nchar(); j++) {
                if (pairs[i][j] != null) {
                    builder.append("<td>").append(pairs[i][j].toString()).append("</td>");
                } else {
                    builder.append("<td></td>"); // Empty cell for null pairs
                }
            }
            builder.append("</tr>\n");
        }

        builder.append("</table>\n");
        builder.append("<p>ntax = ").append(taxa.ntaxa()).append("</p>\n");
        builder.append("</body>\n");
        builder.append("</html>");

        return builder.toString();
    }

    //*** for logging to file ***//

    @Override
    public List<String> getTextForFile() {
        return List.of();
    }

    @Override
    public String getFileType() {
        return ".tsv";
    }
}

