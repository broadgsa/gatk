package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.genotype.BasicVariation;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;


/**
 * 
 * @author aaron 
 * 
 * Class GenotypeUtilsTest
 *
 * a test class for the various function in the genotype utils class
 */
public class GenotypeUtilsTest extends BaseTest {

    private static IndexedFastaSequenceFile seq;
    private static File vcfFile = new File(validationDataLocation + "vcfexample.vcf");

    @BeforeClass
    public static void beforeTests() {
        try {
            seq = new IndexedFastaSequenceFile(new File(oneKGLocation + "reference/human_b36_both.fasta"));
        } catch (FileNotFoundException e) {
            throw new StingException("unable to load the sequence dictionary");
        }
        GenomeLocParser.setupRefContigOrdering(seq);

    }

    /**
     * make sure that the variation is a het
     */
    @Test
    public void isHetTest() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0,1,2);
        BasicVariation var = new BasicVariation("AA", "A", 0, loc,0.0);
        Assert.assertTrue(!GenotypeUtils.isHet(var));
        BasicVariation var2 = new BasicVariation("AG", "A", 0, loc,0.0);
        Assert.assertTrue(GenotypeUtils.isHet(var2));
        BasicVariation var3 = new BasicVariation("GG", "A", 0, loc,0.0);
        Assert.assertTrue(!GenotypeUtils.isHet(var3));
    }

}
