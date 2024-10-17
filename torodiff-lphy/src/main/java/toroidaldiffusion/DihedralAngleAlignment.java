package toroidaldiffusion;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.AugmentedAlignment;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.logger.TextFileFormatted;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * This alignment could contain the internal node sequences.
 * If {@link #internalNodes} is null, then this only contains the tips sequences.
 * If it is not null, then this alignment also contains internal nodes sequences.
 * In this case, the sequence indices for tips range from 0 to (ntaxa - 1),
 * internal node indices follow after the tips, the root node is assigned the final index.
 */
public class DihedralAngleAlignment implements AugmentedAlignment<Pair>, TextFileFormatted {

    // 1st[] is taxa, index is same order as Taxa. When incl. internal nodes, idx is ntaxa * 2 - 1
    // 2nd[] is the site 
    public Pair[][] pairs;
    Taxa taxa;
    List<TimeTreeNode> internalNodes;

    public DihedralAngleAlignment(Taxa taxa, int nchar, List<TimeTreeNode> internalNodes) {
        this.taxa = taxa;
        if (internalNodes == null)
            this.pairs = new Pair[taxa.ntaxa()][nchar];
        else if (internalNodes.size() == taxa.ntaxa() - 1) {
            // include internal nodes
            this.pairs = new Pair[2 * taxa.ntaxa() - 1][nchar];
            this.internalNodes = internalNodes;
        } else
            throw new IllegalArgumentException("The number of internal nodes " + internalNodes.size() +
                    " should be the number of taxa " + taxa.ntaxa() + " - 1 !");
    }

    /**
     * @param sequenceName   check if it is taxon name, if not found, then try the internal node ID,
     *                       if still not found, throw an exception.
     * @param position         site position
     * @return               the state for the given sequence name.
     */
    @Override
    public Pair getState(String sequenceName, int position) {
        int sequenceId = indexOfSequence(sequenceName, taxa, internalNodes);
        if (sequenceId < 0)
            throw new IllegalArgumentException("The sequence name " + sequenceName + " does not exist !");
        return pairs[sequenceId][position];
    }

    @Override
    public Pair getState(int sequenceId, int position) {
        return pairs[sequenceId][position];
    }

    @Override
    public void setState(int sequenceId, int position, Pair state) {
        pairs[sequenceId][position] = state;
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
        // Objects.requireNonNull(taxa).ntaxa()
        String[] taxaNames = taxa.getTaxaNames();
        for (int i = 0; i < Objects.requireNonNull(taxaNames).length; i++) {
            builder.append("  ").append(taxaNames[i]); //getname not tostring
            // taxonIndex may not be same as i
            int taxonIndex = taxa.indexOfTaxon(taxaNames[i]);
            builder.append(" = ").append(Arrays.toString(pairs[taxonIndex]));
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
                // taxonIndex may not be same as i, use getState(taxaNames[i], j)
                if (getState(taxaNames[i], j) != null) {
                    builder.append("<td>").append(getState(taxaNames[i], j).toString()).append("</td>");
                } else {
                    builder.append("<td></td>"); // Empty cell for null pairs
                }
            }
            builder.append("</tr>\n");
        }
        // internal nodes
        for (int i = taxaNames.length; i < pairs.length; i++) {
            builder.append("<tr>");
            builder.append("<td style=\"color:#ADD8E6;\">").append(i).append("</td>");
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
        List<String> lines = new ArrayList<>();

        // Create the header line
        StringBuilder header = new StringBuilder();
        header.append("Taxon");
        for (int j = 0; j < nchar(); j++) {
            header.append("\tsite ").append(j + 1).append("_phi");
            header.append("\tsite ").append(j + 1).append("_psi");
        }
        lines.add(header.toString());

        // Create the data rows for each taxon
        String[] taxaNames = taxa.getTaxaNames();
        for (int i = 0; i < taxaNames.length; i++) {
            StringBuilder row = new StringBuilder();
            row.append(taxaNames[i]);  // Add the taxon name
            for (int j = 0; j < nchar(); j++) {
                if (pairs[i][j] != null) {
                    row.append("\t").append(pairs[i][j].getPhi() != null ? pairs[i][j].getPhi().toString() : "");
                    row.append("\t").append(pairs[i][j].getPsi() != null ? pairs[i][j].getPsi().toString() : "");
                } else {
                    row.append("\t").append("");  // Empty for null phi
                    row.append("\t").append("");  // Empty for null psi
                }
            }
            lines.add(row.toString());
        }

        // internal nodes
        for (int i = taxaNames.length; i < pairs.length; i++) {
            StringBuilder row = new StringBuilder();
            row.append("Internal_Node_").append(i);  // Identifier for internal nodes
            for (int j = 0; j < nchar(); j++) {
                if (pairs[i][j] != null) {
                    row.append("\t").append(pairs[i][j].getPhi() != null ? pairs[i][j].getPhi().toString() : "");
                    row.append("\t").append(pairs[i][j].getPsi() != null ? pairs[i][j].getPsi().toString() : "");
                } else {
                    row.append("\t").append("");  // Empty for null phi
                    row.append("\t").append("");  // Empty for null psi
                }
            }
            lines.add(row.toString());
        }
            //lines.add("");

        return lines;
    }

    @Override
    public String getFileType() {
        return ".tsv";
    }

}

