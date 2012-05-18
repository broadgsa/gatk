/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.report;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

public class GATKReportUnitTest extends BaseTest {
    @Test
    public void testParse() throws Exception {
        String reportPath = validationDataLocation + "exampleGATKReportv2.tbl";
        GATKReport report = new GATKReport(reportPath);
        Assert.assertEquals(report.getVersion(), GATKReportVersion.V1_1);
        Assert.assertEquals(report.getTables().size(), 5);

        GATKReportTableV2 countVariants = report.getTable("CountVariants");
        Assert.assertEquals(countVariants.get(0, "nProcessedLoci"), "63025520");
        Assert.assertEquals(countVariants.get(0, "nNoCalls"), "0");
        Assert.assertEquals(countVariants.get(0, "heterozygosity"), 4.73e-06);

        GATKReportTableV2 validationReport = report.getTable("ValidationReport");
        Assert.assertEquals(validationReport.get(2, "PPV"), Double.NaN);
    }

    @DataProvider(name = "rightAlignValues")
    public Object[][] getRightAlignValues() {
        return new Object[][]{
                new Object[]{null, true},
                new Object[]{"null", true},
                new Object[]{"NA", true},
                new Object[]{"0", true},
                new Object[]{"0.0", true},
                new Object[]{"-0", true},
                new Object[]{"-0.0", true},
                new Object[]{String.valueOf(Long.MAX_VALUE), true},
                new Object[]{String.valueOf(Long.MIN_VALUE), true},
                new Object[]{String.valueOf(Float.MIN_NORMAL), true},
                new Object[]{String.valueOf(Double.MAX_VALUE), true},
                new Object[]{String.valueOf(Double.MIN_VALUE), true},
                new Object[]{String.valueOf(Double.POSITIVE_INFINITY), true},
                new Object[]{String.valueOf(Double.NEGATIVE_INFINITY), true},
                new Object[]{String.valueOf(Double.NaN), true},
                new Object[]{"hello", false}
        };
    }

    @Test(dataProvider = "rightAlignValues")
    public void testIsRightAlign(String value, boolean expected) {
        Assert.assertEquals(GATKReportColumn.isRightAlign(value), expected, "right align of '" + value + "'");
    }

    private GATKReportTableV2 makeBasicTable() {
        GATKReport report = GATKReport.newSimpleReport("TableName", "sample", "value");
        GATKReportTableV2 table = report.getTable("TableName");
        report.addRow("foo.1", "hello");
        report.addRow("foo.2", "world");
        return table;
    }

    @Test
    public void testDottedSampleName() {
        GATKReportTableV2 table = makeBasicTable();
        Assert.assertEquals(table.get(0, "value"), "hello");
        Assert.assertEquals(table.get(1, "value"), "world");
    }

    @Test
    public void testSimpleGATKReport() {
        // Create a new simple GATK report named "TableName" with columns: Roger, is, and Awesome
        GATKReport report = GATKReport.newSimpleReport("TableName", "Roger", "is", "Awesome");

        // Add data to simple GATK report
        report.addRow(12, 23.45, true);
        report.addRow("ans", '3', 24.5);
        report.addRow("hi", "", 2.3);

        // Print the report to console
        //report.print(System.out);

        try {
            File file = createTempFile("GATKReportGatherer-UnitTest", ".tbl");
            //System.out.format("The temporary file" + " has been created: %s%n", file);
            PrintStream ps = new PrintStream(file);
            report.print(ps);
            //System.out.println("File succesfully outputed!");
            GATKReport inputRead = new GATKReport(file);
            //System.out.println("File succesfully read!");
            //inputRead.print(System.out);
            Assert.assertTrue(report.isSameFormat(inputRead));

        } catch (IOException x) {
            System.err.format("IOException: %s%n", x);
        }

    }

    @Test
    public void testGATKReportGatherer() {

        GATKReport report1, report2, report3;
        report1 = new GATKReport();
        report1.addTable("TableName", "Description", 2);
        report1.getTable("TableName").addColumn("colA", "%s");
        report1.getTable("TableName").addColumn("colB", "%c");
        report1.getTable("TableName").set(0, "colA", "NotNum");
        report1.getTable("TableName").set(0, "colB", (char) 64);

        report2 = new GATKReport();
        report2.addTable("TableName", "Description", 2);
        report2.getTable("TableName").addColumn("colA", "%s");
        report2.getTable("TableName").addColumn("colB", "%c");
        report2.getTable("TableName").set(0, "colA", "df3");
        report2.getTable("TableName").set(0, "colB", 'A');

        report3 = new GATKReport();
        report3.addTable("TableName", "Description", 2);
        report3.getTable("TableName").addColumn("colA", "%s");
        report3.getTable("TableName").addColumn("colB", "%c");
        report3.getTable("TableName").set(0, "colA", "df5f");
        report3.getTable("TableName").set(0, "colB", 'c');

        report1.concat(report2);
        report1.concat(report3);

        report1.addTable("Table2", "To contain some more data types", 3);
        GATKReportTableV2 table = report1.getTable("Table2");
        table.addColumn("SomeInt", "%d");
        table.addColumn("SomeFloat", "%.16E");
        table.addColumn("TrueFalse", "%B");
        table.addRowIDMapping("12df", 0);
        table.addRowIDMapping("5f", 1);
        table.addRowIDMapping("RZ", 2);
        table.set("12df", "SomeInt", Byte.MAX_VALUE);
        table.set("12df", "SomeFloat", 34.0);
        table.set("12df", "TrueFalse", true);
        table.set("5f", "SomeInt", Short.MAX_VALUE);
        table.set("5f", "SomeFloat", Double.MAX_VALUE);
        table.set("5f", "TrueFalse", false);
        table.set("RZ", "SomeInt", Long.MAX_VALUE);
        table.set("RZ", "SomeFloat", 535646345.657453464576);
        table.set("RZ", "TrueFalse", true);

        report1.addTable("Table3", "blah", 1, true);
        report1.getTable("Table3").addColumn("a");
        report1.getTable("Table3").addRowIDMapping("q", 2);
        report1.getTable("Table3").addRowIDMapping("5", 3);
        report1.getTable("Table3").addRowIDMapping("573s", 0);
        report1.getTable("Table3").addRowIDMapping("ZZZ", 1);
        report1.getTable("Table3").set("q", "a", "34");
        report1.getTable("Table3").set("5", "a", "c4g34");
        report1.getTable("Table3").set("573s", "a", "fDlwueg");
        report1.getTable("Table3").set("ZZZ", "a", "Dfs");

        try {
            File file = createTempFile("GATKReportGatherer-UnitTest", ".tbl");
            //System.out.format("The temporary file" + " has been created: %s%n", file);
            PrintStream ps = new PrintStream(file);
            report1.print(ps);
            //System.out.println("File succesfully outputed!");
            GATKReport inputRead = new GATKReport(file);
            //System.out.println("File succesfully read!");
            //inputRead.print(System.out);
            Assert.assertTrue(report1.isSameFormat(inputRead));
            Assert.assertTrue(report1.equals(inputRead));

        } catch (IOException x) {
            System.err.format("IOException: %s%n", x);
        }
    }
}
