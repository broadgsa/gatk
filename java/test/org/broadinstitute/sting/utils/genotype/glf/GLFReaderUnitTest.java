package org.broadinstitute.sting.utils.genotype.glf;

import org.broadinstitute.sting.BaseTest;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class GLFReaderUnitTest
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public class GLFReaderUnitTest extends BaseTest {


    // our test file
    static final File glfFile = new File(validationDataLocation + "index_test_likelihoods.glf");
    //static final File glfFile = new File("CALLS.glf");
    static final int finalRecordCount = 484140; // the number of records in the above file
    static final int contigCount = 25;

    /** read in the records from the file */
    @Test
    public void testReadRecords() {
        int recCount = 0;
        List<String> contigs = new ArrayList<String>();
        try {
            GLFReader reader = new GLFReader(glfFile);
            long location = 1;
            while (reader.hasNext()) {
                GLFRecord rec = reader.next();
                if (!contigs.contains(rec.getContig())) {
                    contigs.add(rec.getContig());
                }
                location = rec.getPosition();
                //System.err.println("Record count = " + finalRecordCount + " offset " + rec.offset + " location = " + location + " type = " + rec.getRecordType());
                ++recCount;
            }
        } catch (Exception e) {
            System.err.println("Record count = " + recCount);
            e.printStackTrace();
        }
        Assert.assertEquals(recCount, finalRecordCount);
        Assert.assertEquals(contigs.size(), contigCount);
    }
}
