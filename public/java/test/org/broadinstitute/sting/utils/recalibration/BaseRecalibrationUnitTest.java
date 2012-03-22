package org.broadinstitute.sting.utils.recalibration;

import org.testng.annotations.Test;

import java.io.File;

/**
 * Unit tests for on-the-fly recalibration.
 *
 * @author Mauricio Carneiro
 * @since 3/16/12
 */
public class BaseRecalibrationUnitTest {

    @Test(enabled=true)
    public void testReadingCSV() {
        File csv = new File("public/testdata/exampleCSV.csv");
        BaseRecalibration baseRecalibration = new BaseRecalibration(csv);
        System.out.println("Success");
    }
}
