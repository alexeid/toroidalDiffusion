package toroidaldiffusion.lphybeast.tobeast.generator;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class TimeTreeToBEASTTest {

    final String TREE1 = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";

    @BeforeEach
    void setUp() {

    }

    @Test
    void testNewickLT() {
        TreeParser treeParser = new TreeParser();
        treeParser.initByName("IsLabelledNewick", true, "newick", TREE1);

        assertInternalNodes(treeParser);

    }

    @Test
    void testNewickLT2() {
        TaxonSet taxonset = new TaxonSet();
        Taxon taxon0 = new Taxon("0"); // TODO must start from 1 if there is tree op
        Taxon taxon1 = new Taxon("1");
        Taxon taxon2 = new Taxon("2");
        Taxon taxon3 = new Taxon("3");
        taxonset.initByName("taxon", taxon0, "taxon", taxon1, "taxon", taxon2, "taxon", taxon3);

        TreeParser treeParser = new TreeParser();
        treeParser.initByName("IsLabelledNewick", true,
               "taxonset", taxonset, "newick", TREE1,
                "adjustTipHeights", false); // defeat adjustTipHeights=true

        assertInternalNodes(treeParser);

        //TODO tree op mess up internal node ID ?
        assertEquals(TREE1, treeParser + ";", "Original tree != parsed tree ");
//        System.out.println(treeParser);
    }

    private void assertInternalNodes(TreeParser treeParser) {
        List<Node> internalNodes = treeParser.getInternalNodes();

        for (int i = 0; i < internalNodes.size(); i++) {
            Node node = internalNodes.get(i);

            int nr = node.getNr();
            String id = node.getID();

            assertEquals(Integer.parseInt(id), nr, "i = " + 1);
        }
    }

//    @Test
//    void testNewickLF() {
//        TreeParser treeParser = new TreeParser();
//        treeParser.initByName("IsLabelledNewick", false, "newick", TREE1);
////TODO lack of taxa (Alignment) ?
//
//        assertInternalNodes(treeParser);
//
//    }

}