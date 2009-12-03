package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Dec 2, 2009
 * Time: 1:05:58 AM
 * <p/>
 * some quick tests for the BasicVariation class
 */
public class BasicVariationTest extends BaseTest {
    private static IndexedFastaSequenceFile seq;

    @BeforeClass
    public static void beforeTests() {
        try {
            seq = new IndexedFastaSequenceFile(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        } catch (FileNotFoundException e) {
            throw new StingException("unable to load the sequence dictionary");
        }
        GenomeLocParser.setupRefContigOrdering(seq);

    }

    @Test
    public void testIsBiallelic() {
        BasicVariation var = new BasicVariation("CC", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(!var.isBiallelic());
        BasicVariation var2 = new BasicVariation("CA", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var2.isBiallelic());
        BasicVariation var3 = new BasicVariation("CC", "A", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var3.isBiallelic());
    }

    @Test
    public void testVariantType() {
        // test reference
        BasicVariation var = new BasicVariation("CC", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var.getType() == Variation.VARIANT_TYPE.REFERENCE);

        // test SNP's
        BasicVariation var2 = new BasicVariation("CA", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var2.getType() == Variation.VARIANT_TYPE.SNP);
        BasicVariation var3 = new BasicVariation("AA", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var3.getType() == Variation.VARIANT_TYPE.SNP);

        // test deletions
        BasicVariation var4 = new BasicVariation("", "C", -10, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var4.getType() == Variation.VARIANT_TYPE.DELETION);

        // test insertions
        BasicVariation var5 = new BasicVariation("ACACACACACA", "C", "ACACACACACA".length(), GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var5.getType() == Variation.VARIANT_TYPE.INSERTION);

    }

    @Test(expected = IllegalStateException.class)
    public void testGetAlternativeBaseForSNPNotASNP() {
        // test reference
        BasicVariation var = new BasicVariation("CC", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        var.getAlternativeBaseForSNP();
    }

    @Test(expected = IllegalStateException.class)
    public void testGetAlternativeBaseForSNPFromIndel() {
        // test reference
        BasicVariation var = new BasicVariation("ACACACACACA", "C", "ACACACACACA".length(), GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        var.getAlternativeBaseForSNP();
    }

    @Test(expected = IllegalStateException.class)
    public void testGetAlternativeBaseForSNPFromDel() {
        // test reference
        BasicVariation var = new BasicVariation("", "C", -10, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        var.getAlternativeBaseForSNP();
    }

    @Test
    public void testGetAlternativeBaseForSNP() {
        // test SNP's
        BasicVariation var = new BasicVariation("CA", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertEquals('A', var.getAlternativeBaseForSNP());
        var = new BasicVariation("AC", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertEquals('A', var.getAlternativeBaseForSNP());
        var = new BasicVariation("AA", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertEquals('A', var.getAlternativeBaseForSNP());
    }

    @Test
    public void testGetAlleleList() {
        BasicVariation var = new BasicVariation("CA", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var.getAlleleList().size() == 2);
        Assert.assertTrue(var.getAlleleList().contains("C"));
        Assert.assertTrue(var.getAlleleList().contains("A"));

        var = new BasicVariation("AC", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var.getAlleleList().size() == 2);
        Assert.assertTrue(var.getAlleleList().contains("C"));
        Assert.assertTrue(var.getAlleleList().contains("A"));

        var = new BasicVariation("AA", "C", 0, GenomeLocParser.createGenomeLoc(1, 1, 1), 1.22);
        Assert.assertTrue(var.getAlleleList().size() == 2);
        Assert.assertTrue(var.getAlleleList().get(0).equals("A"));
        Assert.assertTrue(var.getAlleleList().get(1).equals("A"));
    }
}
