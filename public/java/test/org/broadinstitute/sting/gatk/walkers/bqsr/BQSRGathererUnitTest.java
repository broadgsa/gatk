package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.testng.annotations.Test;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Mauricio Carneiro
 * @since 3/7/12
 */
public class BQSRGathererUnitTest {
    RecalibrationArgumentCollection RAC;

    private static File recal1 = new File("public/testdata/exampleCSV.csv");
    private static File recal2 = new File("public/testdata/exampleCSV.2.csv");

    @Test(enabled = false)
    public void testCombineTwoFiles() {
        BQSRGatherer gatherer = new BQSRGatherer();
        List<File> recalFiles = new LinkedList<File> ();
        File output = new File("foo.csv");
        
        recalFiles.add(recal1);
        recalFiles.add(recal2);
        gatherer.gather(recalFiles, output);
    }
}
