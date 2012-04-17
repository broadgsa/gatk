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
        String reportPath = validationDataLocation + "exampleGATKReportv1.tbl";
        GATKReport report = new GATKReport(reportPath);
        Assert.assertEquals(report.getVersion(), GATKReportVersion.V1_0);
        Assert.assertEquals(report.getTables().size(), 5);

        GATKReportTable countVariants = report.getTable("CountVariants");
        Object countVariantsPK = countVariants.getPrimaryKeyByData("CountVariants", "dbsnp", "eval", "none", "all");
        Assert.assertEquals(countVariants.get(countVariantsPK, "nProcessedLoci"), "63025520");
        Assert.assertEquals(countVariants.get(countVariantsPK, "nNoCalls"), "0");
        Assert.assertEquals(countVariants.get(countVariantsPK, "heterozygosity"), 4.73e-06);

        GATKReportTable validationReport = report.getTable("ValidationReport");
        Object validationReportPK = countVariants.getPrimaryKeyByData("CountVariants", "dbsnp", "eval", "none", "novel");
        Assert.assertEquals(validationReport.get(validationReportPK, "PPV"), Double.NaN);
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

    private GATKReportTable makeBasicTable() {
        GATKReport report = GATKReport.newSimpleReport("TableName", "sample", "value");
        GATKReportTable table = report.getTable("TableName");
        report.addRow("foo.1", "hello");
        report.addRow("foo.2", "world");
        return table;
    }

    @Test
    public void testDottedSampleName() {
        GATKReportTable table = makeBasicTable();
        Object pk;

        pk = table.getPrimaryKeyByData("foo.1");
        Assert.assertEquals(table.get(pk, "value"), "hello");

        pk = table.getPrimaryKeyByData("foo.2");
        Assert.assertEquals(table.get(pk, "value"), "world");
    }

    @Test
    public void testFindPrimaryKeyByData() {
        GATKReportTable table = makeBasicTable();
        Assert.assertNotNull(table.findPrimaryKeyByData("foo.1"));
        Assert.assertNotNull(table.findPrimaryKeyByData("foo.1", "hello"));
        Assert.assertNotNull(table.findPrimaryKeyByData("foo.2"));
        Assert.assertNotNull(table.findPrimaryKeyByData("foo.2", "world"));
        Assert.assertNull(table.findPrimaryKeyByData("list", "longer", "than", "column", "count"));
        Assert.assertNull(table.findPrimaryKeyByData("short"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyFindPrimaryKeyByData() {
        GATKReportTable table = makeBasicTable();
        table.findPrimaryKeyByData();
    }

    @Test(expectedExceptions = NullPointerException.class)
    public void testNullFindPrimaryKeyByData() {
        GATKReportTable table = makeBasicTable();
        table.findPrimaryKeyByData((Object[]) null);
    }

    @Test
    public void testSimpleGATKReport() {
        // Create a new simple GATK report named "TableName" with columns: Roger, is, and Awesome
        GATKReport report = GATKReport.newSimpleReport("TableName", "Roger", "is", "Awesome");

        // Add data to simple GATK report
        report.addRow( 12, 23.45, true);
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
        boolean displayPK = false;

        GATKReport report1, report2, report3;
        report1 = new GATKReport();
        report1.addTable("TableName", "Description");
        report1.getTable("TableName").addPrimaryKey("id", displayPK);
        report1.getTable("TableName").addColumn("colA", GATKReportDataType.String.getDefaultValue(), "%s");
        report1.getTable("TableName").addColumn("colB", GATKReportDataType.Character.getDefaultValue(), "%c");
        report1.getTable("TableName").set(1, "colA", "NotNum");
        report1.getTable("TableName").set(1, "colB", (char) 64);

        report2 = new GATKReport();
        report2.addTable("TableName", "Description");
        report2.getTable("TableName").addPrimaryKey("id", displayPK);
        report2.getTable("TableName").addColumn("colA", GATKReportDataType.String.getDefaultValue(), "%s");
        report2.getTable("TableName").addColumn("colB", GATKReportDataType.Character.getDefaultValue(), "%c");
        report2.getTable("TableName").set(2, "colA", "df3");
        report2.getTable("TableName").set(2, "colB", 'A');

        report3 = new GATKReport();
        report3.addTable("TableName", "Description");
        report3.getTable("TableName").addPrimaryKey("id", displayPK);
        report3.getTable("TableName").addColumn("colA", GATKReportDataType.String.getDefaultValue(), "%s");
        report3.getTable("TableName").addColumn("colB", GATKReportDataType.Character.getDefaultValue(), "%c");
        report3.getTable("TableName").set(3, "colA", "df5f");
        report3.getTable("TableName").set(3, "colB", 'c');

        report1.combineWith(report2);
        report1.combineWith(report3);

        report1.addTable("Table2", "To contain some more data types");
        GATKReportTable table = report1.getTable("Table2");
        table.addPrimaryKey("KEY");
        table.addColumn("SomeInt", GATKReportDataType.Integer.getDefaultValue(), true, "%d");
        table.addColumn("SomeFloat", GATKReportDataType.Decimal.getDefaultValue(), true, "%.16E");
        table.addColumn("TrueFalse", false, true, "%B");
        table.set("12df", "SomeInt", Byte.MAX_VALUE);
        table.set("12df", "SomeFloat", 34.0);
        table.set("12df", "TrueFalse", true);
        table.set("5f", "SomeInt", Short.MAX_VALUE);
        table.set("5f", "SomeFloat", Double.MAX_VALUE);
        table.set("5f", "TrueFalse", false);
        table.set("RZ", "SomeInt", Long.MAX_VALUE);
        table.set("RZ", "SomeFloat", 535646345.657453464576);
        table.set("RZ", "TrueFalse", true);

        report1.addTable("Table3", "blah");
        report1.getTable("Table3").addPrimaryKey("HAI");
        report1.getTable("Table3").addColumn("a", true, GATKReportDataType.String.getDefaultFormatString());
        report1.getTable("Table3").set("q", "a", "34");
        report1.getTable("Table3").set("5", "a", "c4g34");
        report1.getTable("Table3").set("573s", "a", "fDlwueg");
        report1.getTable("Table3").set("ZZZ", "a", "Dfs");

        //report1.print(System.out);


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

        //Assert.assertEquals(1,1);

    }
}
