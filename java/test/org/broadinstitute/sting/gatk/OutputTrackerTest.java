package org.broadinstitute.sting.gatk;

import org.junit.Test;
import org.junit.After;
import org.junit.Assert;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner; /**
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

public class OutputTrackerTest {
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
        OutputTracker ot = new OutputTracker(null,null);
        Assert.assertSame("OutputTracker: Output stream incorrectly initialized.", System.out, ot.getOutStream());
        Assert.assertSame("OutputTracker: Error stream incorrectly initialized.", System.err, ot.getErrStream());
    }

    @Test
    public void testOutputStreamAlone() throws FileNotFoundException {
        OutputTracker ot = new OutputTracker(OUTPUT_FILENAME,null);

        final String OUTPUT_TEXT = "out stream test";
        ot.getOutStream().append(OUTPUT_TEXT);

        Scanner outScanner = new Scanner(new File(OUTPUT_FILENAME));
        String outText = outScanner.nextLine();
        Assert.assertFalse("Out stream has too much data", outScanner.hasNext());

        Assert.assertEquals("OutputTracker: Written output is incorrect", outText, OUTPUT_TEXT);
        Assert.assertSame("OutputTracker: Error stream incorrectly initialized.", System.err, ot.getErrStream());
    }

    @Test
    public void testErrorStreamAlone() throws FileNotFoundException {
        OutputTracker ot = new OutputTracker(null,ERROR_FILENAME);

        final String ERROR_TEXT = "err stream test";
        ot.getErrStream().append(ERROR_TEXT);

        Scanner errScanner = new Scanner(new File(ERROR_FILENAME));
        String errText = errScanner.nextLine();
        Assert.assertFalse("Err stream has too much data", errScanner.hasNext());

        Assert.assertSame("OutputTracker: Output stream incorrectly initialized.", System.out, ot.getOutStream());
        Assert.assertEquals("OutputTracker: Written error text is incorrect", errText, ERROR_TEXT);
    }

    @Test
    public void testIndependentStreams() throws FileNotFoundException {
        OutputTracker ot = new OutputTracker(OUTPUT_FILENAME,ERROR_FILENAME);

        final String OUTPUT_TEXT = "out stream test";
        ot.getOutStream().append(OUTPUT_TEXT);

        final String ERROR_TEXT = "err stream test";
        ot.getErrStream().append(ERROR_TEXT);

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
        OutputTracker ot = new OutputTracker(OUTPUT_FILENAME,OUTPUT_FILENAME);
        Assert.assertSame("OutputTracker: FileOutputStreams don't match", ot.getOutFile(), ot.getErrFile());
        Assert.assertSame("OutputTracker: PrintStreams don't match", ot.getOutStream(), ot.getErrStream());
    }
}
