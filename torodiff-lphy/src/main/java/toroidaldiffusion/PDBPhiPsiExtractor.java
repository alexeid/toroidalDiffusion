package toroidaldiffusion;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class PDBPhiPsiExtractor {

    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("Usage: java toroidalDiffusion.PDBPhiPsiExtractor <input PDB file> <output file>");
            return;
        }

        String pdbFilePath = args[0];
        String outputFilePath = args[1];

        PDBFileReader pdbReader = new PDBFileReader();
        try {
            Structure structure = pdbReader.getStructure(pdbFilePath);
            extractPhiPsiAngles(structure, outputFilePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void extractPhiPsiAngles(Structure structure, String outputFilePath) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilePath))) {
            for (Chain chain : structure.getChains()) {
                List<Group> groups = chain.getAtomGroups(GroupType.AMINOACID);
                for (int i = 1; i < groups.size() - 1; i++) {
                    AminoAcid previousAminoAcid = (AminoAcid) groups.get(i - 1);
                    AminoAcid aminoAcid = (AminoAcid) groups.get(i);
                    AminoAcid nextAminoAcid = (AminoAcid) groups.get(i + 1);
                    try {
                        double phi = Calc.getPhi(previousAminoAcid, aminoAcid);
                        double psi = Calc.getPsi(aminoAcid, nextAminoAcid);
                        writer.write(String.format("Residue: %s %s, Phi: %.2f, Psi: %.2f%n",
                                aminoAcid.getPDBName(),
                                aminoAcid.getResidueNumber(),
                                phi,
                                psi));
                    } catch (StructureException e) {
                        // Handle residues where phi/psi can't be calculated
                        writer.write(String.format("Residue: %s %s, Phi: N/A, Psi: N/A%n",
                                aminoAcid.getPDBName(),
                                aminoAcid.getResidueNumber()));
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
