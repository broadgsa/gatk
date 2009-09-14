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
 * 
 * @author aaron 
 * 
 * Class BasicGenotypeTest
 *
 * tests the basic genotype class
 */
public class BasicGenotypeTest extends BaseTest {
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
    public void testBasicGenotypeIsHom() {
        BasicGenotype gt = new BasicGenotype(GenomeLocParser.createGenomeLoc("chr1",1,1),"AA",'A',0);
        Assert.assertTrue(gt.isHom());
        BasicGenotype gt2 = new BasicGenotype(GenomeLocParser.createGenomeLoc("chr1",1,1),"GA",'A',0);
        Assert.assertTrue(!gt2.isHom());
    }

    @Test
    public void testBasicGenotypeIsHet() {
        BasicGenotype gt = new BasicGenotype(GenomeLocParser.createGenomeLoc("chr1",1,1),"AA",'A',0);
        Assert.assertTrue(!gt.isHet());
        BasicGenotype gt2 = new BasicGenotype(GenomeLocParser.createGenomeLoc("chr1",1,1),"GA",'A',0);
        Assert.assertTrue(gt2.isHet());
    }

    @Test
    public void testBasicGenotypeIsVariant() {
        BasicGenotype gt = new BasicGenotype(GenomeLocParser.createGenomeLoc("chr1",1,1),"AA",'A',0);
        Assert.assertTrue(!gt.isVariant('A'));
        BasicGenotype gt2 = new BasicGenotype(GenomeLocParser.createGenomeLoc("chr1",1,1),"GA",'A',0);
        Assert.assertTrue(gt2.isVariant('A'));
        BasicGenotype gt3 = new BasicGenotype(GenomeLocParser.createGenomeLoc("chr1",1,1),"TT",'A',0);
        Assert.assertTrue(gt3.isVariant('A'));
    }
}
