package org.broadinstitute.sting.gatk.refdata;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.genotype.glf.SinglePointCall;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: aaronmckenna
 * Date: Jul 15, 2009
 * Time: 12:18:50 AM
 */
public class RodGLFTest extends BaseTest {
    static final File glfFile = new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/index_test_likelihoods.glf");
    static final int finalRecordCount = 484140; // the number of records in the above file
    static final int contigCount = 25;
    static final String ref = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";

    @BeforeClass
    public static void before() {
        ReferenceSequenceFile r = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(ref));
        GenomeLocParser.setupRefContigOrdering(r);
    }

    static class TestRODClass implements GLFRODIterator {
        private int count = 0;
        private int readCount = 0;

        public TestRODClass(int count) {
            this.count = count;
        }

        @Override
        public boolean hasNext() {
            if (count <= 0) {
                return false;
            }
            return true;
        }

        @Override
        public RodGLF next() {
            RodGLF glf = new RodGLF("Test");
            glf.mRecord = new SinglePointCall('A', readCount, 0, (short) 0, new double[]{0.0, 0.0});
            --count;
            readCount++;
            return glf;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("BOO");
        }
    }

    //@Test
    public void testNext() {
        RodGLF.setGLFWrapper(new RodGLFTest.TestRODClass(10));
        Iterator<RodGLF> iter = RodGLF.createIterator("test", new File(""));
        int counter = 0;
        while (iter.hasNext()) {
            iter.next();
            counter++;
        }
        assertEquals(10, counter);
    }

    @Test
    public void testRodCount() {

        Iterator<RodGLF> iter = RodGLF.createIterator("test", glfFile);
        int counter = 0;
        while (iter.hasNext()) {
            RodGLF glf = iter.next();
            counter++;
        }
        assertEquals(finalRecordCount, counter);
    }

    @Test
    public void testRodLocations() {

        Iterator<RodGLF> iter = RodGLF.createIterator("test", glfFile);
        GenomeLoc loc = null;
        while (iter.hasNext()) {
            RodGLF glf = iter.next();
            if (loc != null) {
                if (glf.getLocation().isBefore(loc)) {
                    fail("locations in the GLF came out of order loc = " + loc.toString() + " new loc = " + glf.getLocation().toString());
                }
            }
            loc = glf.getLocation();
        }        
    }


}
