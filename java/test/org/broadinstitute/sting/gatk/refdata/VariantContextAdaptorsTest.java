package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.glf.GLFSingleCall;
import org.broadinstitute.sting.utils.genotype.glf.GLFWriter;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class VariantContextAdaptorsTest
 *
 * This test class exists to test input -> output of variant formats
 * run through the VariantContext class.  
 */
public class VariantContextAdaptorsTest extends BaseTest {
    public static IndexedFastaSequenceFile seq = null;

    @BeforeClass
    public static void beforeClass() {
        try {
            seq = new IndexedFastaSequenceFile(new File(oneKGLocation + "/reference/human_b36_both.fasta")); // TODO: make human reference use BaseTest
        } catch (FileNotFoundException e) {
            Assert.fail("Unable to load reference " + oneKGLocation + "/reference/human_b36_both.fasta");
        }
        GenomeLocParser.setupRefContigOrdering(seq);
    }


    /**
     * this test takes a known GLF file, reads in the records (storing them into an array),
     * and creates VariantContext records.  These VC records are then outputted through a genotype writer,
     * and then read back in off of disk and compared to the original records.  This way we are positive all
     * the information that encodes a GLF makes it into the VC and then out to disk.
     */
    @Test
    public void testVariantContextGLFToGLF() {

        // our input and output files
        File referenceFile = new File(validationDataLocation + "/well_formed.glf");       // our known good GLF
        File tempFile = new File("temp.glf");                   // our temporary GLF output -> input file
        tempFile.deleteOnExit();                                // delete when we're done

        // create our genotype writer for GLFs
        GenotypeWriter gw = GenotypeWriterFactory.create(GenotypeWriterFactory.GENOTYPE_FORMAT.GLF,tempFile);
        ((GLFWriter)gw).writeHeader("");

        RodGLF glf = new RodGLF("myROD"); // now cycle the input file to the output file
        try {
            glf.initialize(referenceFile);
        } catch (FileNotFoundException e) {
            Assert.fail("Unable to open GLF file" + referenceFile);
        }

        // buffer the records we see
        List<GLFSingleCall> records = new ArrayList<GLFSingleCall>();

        // while we have records, make a Variant Context and output it to a GLF file
        while (glf.hasNext()) {
            glf.next();
            records.add((GLFSingleCall)glf.mRecord); // we know they're all single calls in the reference file
            VariantContext vc = VariantContextAdaptors.toVariantContext("GLF",glf);
            gw.addCall(vc);
        }
        gw.close(); // close the file


        // now reopen the file with the temp GLF file and read it back in, compare against what we first stored
        glf = new RodGLF("myROD");
        try {
            glf.initialize(tempFile);
        } catch (FileNotFoundException e) {
            Assert.fail("Unable to open GLF file" + tempFile);
        }
        
        // buffer the new records we see
        List<GLFSingleCall> records2 = new ArrayList<GLFSingleCall>();

        // while we have records, make a Variant Context and output it to a GLF file
        while (glf.hasNext()) {
            glf.next();
            records2.add((GLFSingleCall)glf.mRecord); // we know they're all single calls in the reference file
        }

        // compare sizes
        Assert.assertEquals("The input GLF file doesn't contain the same number of records as we saw in the first file", records.size(),records2.size());

        // now compare each record
        for (int x = 0; x < records.size(); x++)
            Assert.assertTrue("GLF Records were not preserved when cycling them to and from disc", records.get(x).equals(records2.get(x)));
    }
}
