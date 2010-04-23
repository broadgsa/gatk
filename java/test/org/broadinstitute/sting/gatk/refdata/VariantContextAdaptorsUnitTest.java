package org.broadinstitute.sting.gatk.refdata;

import edu.mit.broad.picard.genotype.geli.GeliFileReader;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.util.AsciiLineReader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.geli.GeliGenotypeWriter;
import org.broadinstitute.sting.utils.genotype.geli.GeliTextWriter;
import org.broadinstitute.sting.utils.genotype.glf.GLFSingleCall;
import org.broadinstitute.sting.utils.genotype.glf.GLFWriter;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class VariantContextAdaptorsUnitTest
 *
 * This test class exists to test input -> output of variant formats
 * run through the VariantContext class.  
 */
public class VariantContextAdaptorsUnitTest extends BaseTest {
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
            gw.addCall(vc,null);
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

    /**
     * this test takes a known Geli file, reads in the records (storing them into an array),
     * and creates VariantContext records.  These VC records are then outputted through a genotype writer,
     * and then read back in off of disk and compared to the original records.  This way we are positive all
     * the information that encodes a Geli makes it into the VC and then out to disk.
     */
    @Test
    public void testVariantContextGeliTextToGeliText() {

        // our input and output files
        File knownFile = new File(validationDataLocation + "/well_formed.geli");        // our known good geli
        File tempFile = new File("temp.geli");                                          // our temporary geli output -> input file
        tempFile.deleteOnExit();                                                        // delete when we're done

        // create our genotype writer for GLFs
        GenotypeWriter gw = GenotypeWriterFactory.create(GenotypeWriterFactory.GENOTYPE_FORMAT.GELI,tempFile);
        ((GeliTextWriter)gw).writeHeader(null);  // the write header command ignores the parameter

        RodGeliText geliText = new RodGeliText("myROD"); // now cycle the input file to the output file        

        // buffer the records we see
        List<RodGeliText> records = new ArrayList<RodGeliText>();

        // a little more complicated than the GLF example, we have to read the file in
        AsciiLineReader reader = createReader(knownFile);

        // get the first real line (non-header)
        String line = cleanHeaderFromFile(reader);

        // while we have records, make a Variant Context and output it to a GLF file
        while (line != null && line != "") {
            parseGeliText(geliText, line);
            records.add(geliText); // we know they're all single calls in the reference file
            VariantContext vc = VariantContextAdaptors.toVariantContext("Geli",geliText);
            gw.addCall(vc,null);
            line = readLine(reader);
        }
        gw.close(); // close the file
        reader.close();

        // now reopen the file with the temp GLF file and read it back in, compare against what we first stored
        geliText = new RodGeliText("myROD");

        // buffer the new records we see
        List<RodGeliText> records2 = new ArrayList<RodGeliText>();

        reader = createReader(tempFile);
        // get the first real line (non-header)
        line = cleanHeaderFromFile(reader);

        // while we have records, make a Variant Context and output it to a GLF file
        while (line != null && line != "") {
            parseGeliText(geliText,line);
            records2.add(geliText); // we know they're all single calls in the reference file
            line = readLine(reader);
            
        }
        gw.close(); // close the file
        reader.close();
        // compare sizes
        Assert.assertEquals("The input GLF file doesn't contain the same number of records as we saw in the first file", records.size(),records2.size());
        // now compare each record
        for (int x = 0; x < records.size(); x++)
            Assert.assertTrue("GLF Records were not preserved when cycling them to and from disc", records.get(x).equals(records2.get(x)));
    }

    /**
     * this test takes a known GeliBinary file, reads in the records (storing them into an array),
     * and creates VariantContext records.  These VC records are then outputted through a genotype writer,
     * and then read back in off of disk and compared to the original records.  This way we are positive all
     * the information that encodes a Geli makes it into the VC and then out to disk.
     */
    @Test
    public void testVariantContextGeliBinaryToGeliBinary() {

        // our input and output files
        File knownFile = new File(validationDataLocation + "large_test_geli_binary.geli");         // our known good geli
        File tempFile = new File("temp_binary.geli");                                                      // our temporary geli output -> input file
        //tempFile.deleteOnExit();                                                                    // delete when we're done

        // create our genotype writer for GLFs
        GenotypeWriter gw = GenotypeWriterFactory.create(GenotypeWriterFactory.GENOTYPE_FORMAT.GELI_BINARY,tempFile);

        // buffer the records we see
        List<rodGELI> records = new ArrayList<rodGELI>();

        // a little more complicated than the GLF example, we have to read the file in
        GeliFileReader reader = new GeliFileReader(knownFile);

        ((GeliGenotypeWriter)gw).writeHeader(reader.getFileHeader());  // the write header command ignores the parameter

        CloseableIterator<GenotypeLikelihoods> iterator = reader.iterator();
        // while we have records, make a Variant Context and output it to a GeliBinary file
        while (iterator.hasNext()) {
            rodGELI gel = new rodGELI("myROD",iterator.next());
            records.add(gel);
            VariantContext vc = VariantContextAdaptors.toVariantContext("myROD",gel);
            if (vc != null) gw.addCall(vc,null);
        }
        iterator.close();
        gw.close(); // close the file
        reader.close();

        // buffer the new records we see
        List<rodGELI> records2 = new ArrayList<rodGELI>();

        reader = new GeliFileReader(tempFile);
        iterator = reader.iterator();
        while (iterator.hasNext()) {
            rodGELI gel = new rodGELI("myROD",iterator.next());
            records2.add(gel);
        }
        iterator.close();
        reader.close();
        // compare sizes
        Assert.assertEquals("The input GeliBinary file doesn't contain the same number of records as we saw in the first file", records.size(),records2.size());
        // now compare each record
        for (int x = 0; x < records.size(); x++) {
            Assert.assertTrue("GeliBinary Records were not preserved when cycling them to and from disc", records.get(x).getGeliRecord().equals(records2.get(x).getGeliRecord()));
        }
    }


    /**
     * parse the geli given a line representation
     * @param geliText the object to parse into
     * @param line the line to parse with
     */
    private void parseGeliText(RodGeliText geliText, String line) {
        boolean parsed = false;
        try {
            parsed = geliText.parseLine(null,line.split(TabularROD.DEFAULT_DELIMITER_REGEX));
        } catch (IOException e) {
            Assert.fail("IOException: " + e.getMessage());
        }
        if (!parsed) Assert.fail("Unable to parse line" + line);
    }

    /**
     * clean header lines form a geli file
     * @param reader the reader
     * @return the first line that's not a header
     */
    private String cleanHeaderFromFile(AsciiLineReader reader) {
        String line = "#";
        while (line != null && line.startsWith("#"))
            line = readLine(reader);
        return line;
    }

    /**
     * create a reader, given a file and the reader object
     * @param knownFile the file to read in
     * @return AsciiLineReader the ascii line reader
     */
    private AsciiLineReader createReader(File knownFile) {
        AsciiLineReader reader = null;
        try {
            reader = new AsciiLineReader(new FileInputStream(knownFile));
        } catch (FileNotFoundException e) {
            Assert.fail("File not found: " + knownFile);
        }
        return reader;
    }

    /**
     * read a line from the specified reader.  A method to make the above code look cleaner
     * @param reader the ascii line reader
     * @return a line of text
     */
    public String readLine(AsciiLineReader reader) {
        try {
            String line = reader.readLine();
            return line;
        } catch (IOException e) {
            return null;
        }
    }
}
