package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTableV2;
import org.testng.Assert;
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

    private static File recal = new File("public/testdata/exampleGRP.grp");

    //todo -- this test doesnt work because the primary keys in different tables are not the same. Need to either implement "sort" for testing purposes on GATKReport or have a sophisticated comparison measure
    @Test(enabled = false)
    public void testCombineSimilarFiles() {
        BQSRGatherer gatherer = new BQSRGatherer();
        List<File> recalFiles = new LinkedList<File> ();
        File output = new File("foo.grp");
        recalFiles.add(recal);
        recalFiles.add(recal);
        gatherer.gather(recalFiles, output);

        GATKReport originalReport = new GATKReport(recal);
        GATKReport calculatedReport = new GATKReport(output);        
        for (GATKReportTableV2 originalTable : originalReport.getTables()) {
            GATKReportTableV2 calculatedTable = calculatedReport.getTable(originalTable.getTableName());
            List<String> columnsToTest = new LinkedList<String>();
            columnsToTest.add(RecalDataManager.NUMBER_OBSERVATIONS_COLUMN_NAME);
            columnsToTest.add(RecalDataManager.NUMBER_ERRORS_COLUMN_NAME);
            if (originalTable.getTableName().equals(RecalDataManager.ARGUMENT_REPORT_TABLE_TITLE)) {                    // these tables must be IDENTICAL
                columnsToTest.add(RecalDataManager.ARGUMENT_VALUE_COLUMN_NAME);
                testTablesWithColumnsAndFactor(originalTable, calculatedTable, columnsToTest, 1);
            }
            
            else if (originalTable.getTableName().equals(RecalDataManager.QUANTIZED_REPORT_TABLE_TITLE)) {
                columnsToTest.add(RecalDataManager.QUANTIZED_COUNT_COLUMN_NAME);
                testTablesWithColumnsAndFactor(originalTable, calculatedTable, columnsToTest, 2);
            }
            
            else if (originalTable.getTableName().startsWith("RecalTable")) {
                testTablesWithColumnsAndFactor(originalTable, calculatedTable, columnsToTest, 2);
            }                                              
        }
    }

    /**
     * Common testing functionality given the columns to test and the multiplication factor to the expected result
     *
     * @param original the original table
     * @param calculated the calculated table
     * @param columnsToTest list of columns to test. All columns will be tested with the same criteria (equality given factor)
     * @param factor 1 to test for equality, any other value to multiply the original value and match with the calculated
     */
    private void testTablesWithColumnsAndFactor(GATKReportTableV2 original, GATKReportTableV2 calculated, List<String> columnsToTest, int factor) {
        for (int row = 0; row < original.getNumRows(); row++ ) {
            for (String column : columnsToTest) {
                Object actual = calculated.get(new Integer(row), column);
                Object expected = original.get(row, column);

                if (factor != 1) {
                    if (expected instanceof Double)
                        expected = (Double) expected * factor;
                    else if (expected instanceof Long)
                        expected = (Long) expected * factor;
                    else if (expected instanceof Integer)
                        expected = (Integer) expected * factor;
                    else if (expected instanceof Byte) {
                        expected = (Byte) expected * factor;
                    }
                }
                Assert.assertEquals(actual, expected, "Row: " + row + " Original Table: " + original.getTableName() + " Calc Table: " + calculated.getTableName());
            }
        }
        
    }
}
