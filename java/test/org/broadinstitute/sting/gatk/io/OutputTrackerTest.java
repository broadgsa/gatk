package org.broadinstitute.sting.gatk.io;

import org.junit.Test;
import org.junit.After;
import org.junit.Assert;
import org.broadinstitute.sting.utils.io.RedirectingOutputStream;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.io.DirectOutputTracker;
import org.broadinstitute.sting.gatk.io.stubs.OutputStreamStub;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner; 
/**
 * User: hanna
 * Date: Apr 30, 2009
 * Time: 10:20:18 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * For the file opening and closing mechanisms.
 */

public class OutputTrackerTest extends BaseTest {
    public static final String OUTPUT_FILENAME = "OutputTrackerTest.out";
    public static final String ERROR_FILENAME = "OutputTrackerTest.err";

    @After
    public void cleanupTestFiles() {
        File outFile = new File( OUTPUT_FILENAME );
        File errFile = new File( ERROR_FILENAME );

        if( outFile.exists() ) outFile.delete();
        if( errFile.exists() ) errFile.delete();
    }

    @Test
    public void testNullInputs() {
        OutputTracker ot = new DirectOutputTracker();
        ot.initializeCoreIO(null,null);

        Assert.assertNotNull("OutputTracker: Output stream is null.", ot.outStub );
        Assert.assertNotNull("OutputTracker: Error stream is null.", ot.errStub );

        OutputStreamStub outStream = ot.outStub;
        Assert.assertNull("OutputTracker: Output file incorrectly initialized.", outStream.getOutputFile());
        Assert.assertSame("OutputTracker: Output stream incorrectly initialized.", System.out, outStream.getOutputStream());

        OutputStreamStub errStream = ot.errStub;        
        Assert.assertNull("OutputTracker: Error file incorrectly initialized.", errStream.getOutputFile());
        Assert.assertSame("OutputTracker: Error stream incorrectly initialized.", System.err, errStream.getOutputStream());
    }

    @Test
    public void testOutputStreamAlone() throws FileNotFoundException {
        OutputTracker ot = new DirectOutputTracker();
        ot.initializeCoreIO(OUTPUT_FILENAME,null);

        final String OUTPUT_TEXT = "out stream test";
        PrintWriter outWriter = new PrintWriter(ot.outStub);
        outWriter.append(OUTPUT_TEXT);
        outWriter.close();

        Scanner outScanner = new Scanner(new File(OUTPUT_FILENAME));
        String outText = outScanner.nextLine();
        Assert.assertFalse("Out stream has too much data", outScanner.hasNext());

        Assert.assertEquals("OutputTracker: Written output is incorrect", outText, OUTPUT_TEXT);

        OutputStreamStub errStream = ot.errStub;
        Assert.assertNull("OutputTracker: Error file incorrectly initialized.", errStream.getOutputFile());
        Assert.assertSame("OutputTracker: Error stream incorrectly initialized.", System.err, errStream.getOutputStream());
    }

    //@Test
    public void testErrorStreamAlone() throws FileNotFoundException {
        OutputTracker ot = new DirectOutputTracker();
        ot.initializeCoreIO(null,ERROR_FILENAME);

        OutputStreamStub outStream = ot.outStub;
        Assert.assertNull("OutputTracker: Output file incorrectly initialized.", outStream.getOutputFile());
        Assert.assertSame("OutputTracker: Output stream incorrectly initialized.", System.out, outStream.getOutputStream());

        final String ERROR_TEXT = "err stream test";
        PrintWriter errWriter = new PrintWriter(ot.errStub);
        errWriter.append(ERROR_TEXT);
        errWriter.close();

        Scanner errScanner = new Scanner(new File(ERROR_FILENAME));
        String errText = errScanner.nextLine();
        Assert.assertFalse("Err stream has too much data", errScanner.hasNext());

        Assert.assertEquals("OutputTracker: Written error text is incorrect", errText, ERROR_TEXT);
    }

    //@Test
    public void testIndependentStreams() throws FileNotFoundException {
        OutputTracker ot = new DirectOutputTracker();
        ot.initializeCoreIO(OUTPUT_FILENAME,ERROR_FILENAME);

        final String OUTPUT_TEXT = "out stream test";
        PrintWriter outWriter = new PrintWriter(ot.outStub);
        outWriter.append(OUTPUT_TEXT);
        outWriter.close();

        final String ERROR_TEXT = "err stream test";
        PrintWriter errWriter = new PrintWriter(ot.errStub);
        errWriter.append(ERROR_TEXT);
        errWriter.close();

        Scanner outScanner = new Scanner(new File(OUTPUT_FILENAME));
        String outText = outScanner.nextLine();
        Assert.assertFalse("Out stream has too much data", outScanner.hasNext());

        Scanner errScanner = new Scanner(new File(ERROR_FILENAME));
        String errText = errScanner.nextLine();
        Assert.assertFalse("Err stream has too much data", errScanner.hasNext());

        Assert.assertEquals("OutputTracker: Written output text is incorrect", outText, OUTPUT_TEXT);
        Assert.assertEquals("OutputTracker: Written error text is incorrect", errText, ERROR_TEXT);
    }

    @Test
    public void testIdenticalInputsGetIdenticalResults() {
        OutputTracker ot = new DirectOutputTracker();
        ot.initializeCoreIO(OUTPUT_FILENAME,OUTPUT_FILENAME);

        Assert.assertNotNull("OutputTracker: Output stream is null.", ot.outStub );
        Assert.assertNotNull("OutputTracker: Error stream is null.", ot.errStub );

        OutputStreamStub outStream = ot.outStub;
        OutputStreamStub errStream = ot.errStub;

        Assert.assertSame("OutputTracker: files don't match", outStream.getOutputFile(), errStream.getOutputFile());
        Assert.assertSame("OutputTracker: streams don't match", outStream.getOutputStream(), errStream.getOutputStream());
    }
}
