package org.broadinstitute.sting.gatk.refdata;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.junit.Assert;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Jul 15, 2009
 * Time: 12:18:50 AM
 * <p/>
 * These tests work upon a very small data set, with the following samtools glfview dump:
 * <p/>
 * chrM	1	A	5	20	0	0	127	127	127	254	254	254	254	254	254
 * chrM	2	A	5	20	0	254	254	254	127	254	254	127	254	127	0
 * chrM	3	A	5	20	0	254	127	254	254	0	127	127	254	254	254
 * <p/>
 * You'll notice that the first is a hom ref, and the other two are hom alt SNP's
 */
public class RodGLFTest extends BaseTest {
    static final File glfFile = new File(validationDataLocation + "glfTestFile.glf");
    static final int finalRecordCount = 100; // the number of records in the above file
    static final int contigCount = 1;
    static final String ref = oneKGLocation + "reference/human_b36_both.fasta";
    static ReferenceSequenceFile r;
    private RodGLF iter = null;

    @BeforeClass
    public static void before() {
        r = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(ref));
        GenomeLocParser.setupRefContigOrdering(r);
    }

    @Before
    public void setup() {
        iter = RodGLF.createIterator("test", glfFile);
    }

    @Test
    public void testRodCount() {
        int counter = 0;
        while (iter.hasNext()) {
            RodGLF glf = iter.next();
            counter++;
        }
        assertEquals(finalRecordCount, counter);
    }



    @Test
    public void testIsSNP() {
        RodGLF glf = iter.next();
        Assert.assertFalse(iter.isSNP());
        glf = iter.next();
        Assert.assertTrue(iter.isSNP());
        glf = iter.next();
        Assert.assertTrue(iter.isSNP());
    }

    @Test
    public void testIsReference() {
        RodGLF glf = iter.next();
        Assert.assertTrue(iter.isReference());
        glf = iter.next();
        Assert.assertFalse(iter.isReference());
        glf = iter.next();
        Assert.assertFalse(iter.isReference());
    }
    
    @Test(expected = IllegalStateException.class)
    public void testGetAltSnpFWDIllegalException() {
        RodGLF glf = iter.next();
        iter.getAlternativeBaseForSNP();
    }




    @Test
    public void testCompareTo() {        
        RodGLF iter2 = RodGLF.createIterator("test", glfFile);
        RodGLF glf = iter.next();
        RodGLF glf2 = iter2.next();
        assertEquals(0, glf.compareTo(glf2));
        glf2 = iter2.next();
        assertEquals(-1, glf.compareTo(glf2));
        assertEquals(1, glf2.compareTo(glf));

    }

    @Test
    public void testGetAltSnpFWD() {
        RodGLF glf = iter.next();
        glf = iter.next();
        Assert.assertEquals('C', iter.getAlternativeBaseForSNP());
    }

    @Test
    public void testGetRefSnpFWD() {
        RodGLF glf = iter.next();
        glf = iter.next();
        glf = iter.next();
        Assert.assertEquals('A', iter.getReferenceForSNP());
    }

   

    /**
     * move to the second and third bases, and check that the
     * alternate bases are correct.
     *
    @Test
    public void testGetAltBasesFWD() {
        RodGLF glf = iter.next();
        glf = iter.next();
        Assert.assertTrue("GT".equals(iter.getAltBasesFWD()));
        glf = iter.next();
        Assert.assertTrue("CT".equals(iter.getAltBasesFWD()));

    }

    @Test
    public void testRodLocations() {
        GenomeLoc loc = null;
        while (iter.hasNext()) {
            RodGLF glf = iter.next();
            if (loc != null) {
                if (iter.getLocation().isBefore(loc)) {
                    Assert.fail("locations in the GLF came out of order loc = " + loc.toString() + " new loc = " + iter.getLocation().toString());
                }
            }
            loc = iter.getLocation();
        }
    }

    //@Test
    /**
     * create the example glf file for the test, you can uncomment the above test line to have this
     * test run, regenerating the file.
     *
    public void createRodFile() {
        GenotypeWriter writer = new GLFWriter("", new File("glfTestFile.glf"));
        int location = 1;
        int x = 0;
        writer.addGenotypeCall(r.getSequenceDictionary().getSequence(0), 1, 20, 'A', 5, createLikelihood('A'));
        writer.addGenotypeCall(r.getSequenceDictionary().getSequence(0), 2, 20, 'A', 5, createLikelihood('T'));
        writer.addGenotypeCall(r.getSequenceDictionary().getSequence(0), 3, 20, 'A', 5, createLikelihood('C'));
        writer.close();
    }*/

    /**
     * create a likelihood object, given the appropriate reference base
     *
     * @param ref the reference base
     *
     * @return the likelihood object
     *
    private LikelihoodObject createLikelihood(char ref) {
        ArrayList<Double> vals = new ArrayList<Double>();
        for (LikelihoodObject.GENOTYPE type : LikelihoodObject.GENOTYPE.values()) {
            double x = (type.toString().charAt(0) == ref) ? 0 : 127 - (10 * Math.random());
            x += (type.toString().charAt(1) == ref) ? 0 : 127 - (10 * Math.random());
            vals.add(x);
        }
        double ret[] = new double[vals.size()];
        for (int x = 0; x < vals.size(); x++) {
            ret[x] = vals.get(x);
        }
        return new LikelihoodObject(ret, LikelihoodObject.LIKELIHOOD_TYPE.NEGATIVE_LOG);
    }*/


    /**
     * just make sure that we do get a string back, and no exceptions are thrown
     */
    @Test
    public void testToString() {
        RodGLF glf = iter.next();
        iter.toString();
    }
}