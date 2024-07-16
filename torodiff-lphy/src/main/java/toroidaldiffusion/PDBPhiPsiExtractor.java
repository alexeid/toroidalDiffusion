package toroidaldiffusion;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.secstruc.SecStrucInfo;
import org.biojava.nbio.structure.secstruc.SecStrucTools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class PDBPhiPsiExtractor {

    private static final Map<String, Double> HYDROPHOBICITY_MAP = new HashMap<>();
    private static final Map<String, String> THREE_TO_ONE_MAP = new HashMap<>();

    static {
        // Kyte-Doolittle hydrophobicity scale
        HYDROPHOBICITY_MAP.put("A", 1.8);
        HYDROPHOBICITY_MAP.put("C", 2.5);
        HYDROPHOBICITY_MAP.put("D", -3.5);
        HYDROPHOBICITY_MAP.put("E", -3.5);
        HYDROPHOBICITY_MAP.put("F", 2.8);
        HYDROPHOBICITY_MAP.put("G", -0.4);
        HYDROPHOBICITY_MAP.put("H", -3.2);
        HYDROPHOBICITY_MAP.put("I", 4.5);
        HYDROPHOBICITY_MAP.put("K", -3.9);
        HYDROPHOBICITY_MAP.put("L", 3.8);
        HYDROPHOBICITY_MAP.put("M", 1.9);
        HYDROPHOBICITY_MAP.put("N", -3.5);
        HYDROPHOBICITY_MAP.put("P", -1.6);
        HYDROPHOBICITY_MAP.put("Q", -3.5);
        HYDROPHOBICITY_MAP.put("R", -4.5);
        HYDROPHOBICITY_MAP.put("S", -0.8);
        HYDROPHOBICITY_MAP.put("T", -0.7);
        HYDROPHOBICITY_MAP.put("V", 4.2);
        HYDROPHOBICITY_MAP.put("W", -0.9);
        HYDROPHOBICITY_MAP.put("Y", -1.3);

        // Three-letter to one-letter amino acid codes
        THREE_TO_ONE_MAP.put("ALA", "A");
        THREE_TO_ONE_MAP.put("CYS", "C");
        THREE_TO_ONE_MAP.put("ASP", "D");
        THREE_TO_ONE_MAP.put("GLU", "E");
        THREE_TO_ONE_MAP.put("PHE", "F");
        THREE_TO_ONE_MAP.put("GLY", "G");
        THREE_TO_ONE_MAP.put("HIS", "H");
        THREE_TO_ONE_MAP.put("ILE", "I");
        THREE_TO_ONE_MAP.put("LYS", "K");
        THREE_TO_ONE_MAP.put("LEU", "L");
        THREE_TO_ONE_MAP.put("MET", "M");
        THREE_TO_ONE_MAP.put("ASN", "N");
        THREE_TO_ONE_MAP.put("PRO", "P");
        THREE_TO_ONE_MAP.put("GLN", "Q");
        THREE_TO_ONE_MAP.put("ARG", "R");
        THREE_TO_ONE_MAP.put("SER", "S");
        THREE_TO_ONE_MAP.put("THR", "T");
        THREE_TO_ONE_MAP.put("VAL", "V");
        THREE_TO_ONE_MAP.put("TRP", "W");
        THREE_TO_ONE_MAP.put("TYR", "Y");
    }

    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("Usage: java toroidalDiffusion.PDBPhiPsiExtractor <input PDB file> <output file>");
            return;
        }

        String pdbFilePath = args[0];
        String outputFilePath = args[1];

        FileParsingParameters params = new FileParsingParameters();
        params.setParseSecStruc(true);
        AtomCache cache = new AtomCache();
        cache.setFileParsingParams(params);

        try {
            Structure structure = cache.getStructure(pdbFilePath);
            extractPhiPsiAngles(structure, outputFilePath);
        } catch (IOException | StructureException e) {
            e.printStackTrace();
        }
    }

    public static void extractPhiPsiAngles(Structure structure, String outputFilePath) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilePath))) {
            // Write CSV header
            writer.write("Chain,Residue,ResidueNumber,Phi,Psi,Hydrophobicity,SecondaryStructure\n");

            Map<Group, String> secondaryStructureMap = getSecondaryStructureMap(structure);

            for (Chain chain : structure.getChains()) {
                List<Group> groups = chain.getAtomGroups(GroupType.AMINOACID);
                for (int i = 1; i < groups.size() - 1; i++) {
                    AminoAcid previousAminoAcid = (AminoAcid) groups.get(i - 1);
                    AminoAcid aminoAcid = (AminoAcid) groups.get(i);
                    AminoAcid nextAminoAcid = (AminoAcid) groups.get(i + 1);
                    String residueName = aminoAcid.getPDBName();
                    String oneLetterCode = THREE_TO_ONE_MAP.getOrDefault(residueName, "");
                    double hydrophobicity = HYDROPHOBICITY_MAP.getOrDefault(oneLetterCode, 0.0);
                    String secondaryStructure = secondaryStructureMap.getOrDefault(aminoAcid, "coil");
                    try {
                        double phi = Calc.getPhi(previousAminoAcid, aminoAcid);
                        double psi = Calc.getPsi(aminoAcid, nextAminoAcid);
                        writer.write(String.format("%s,%s,%s,%.2f,%.2f,%.2f,%s%n",
                                chain.getId(),
                                residueName,
                                aminoAcid.getResidueNumber().toString(),
                                phi,
                                psi,
                                hydrophobicity,
                                secondaryStructure));
                    } catch (StructureException e) {
                        // Handle residues where phi/psi can't be calculated
                        writer.write(String.format("%s,%s,%s,N/A,N/A,%.2f,%s%n",
                                chain.getId(),
                                residueName,
                                aminoAcid.getResidueNumber().toString(),
                                hydrophobicity,
                                secondaryStructure));
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Map<Group, String> getSecondaryStructureMap(Structure structure) {
        Map<Group, String> secondaryStructureMap = new HashMap<>();
        List<SecStrucInfo> secStrucInfos = SecStrucTools.getSecStrucInfo(structure);
        for (SecStrucInfo secStrucInfo : secStrucInfos) {
            Group group = secStrucInfo.getGroup();
            secondaryStructureMap.put(group, secStrucInfo.getType().toString());
        }
        return secondaryStructureMap;
    }
}
