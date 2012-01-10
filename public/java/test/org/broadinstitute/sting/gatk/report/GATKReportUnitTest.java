/*
 * Copyright (c) 2011, The Broad Institute
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

public class GATKReportUnitTest extends BaseTest {
    @Test
    public void testParse() throws Exception {
        String reportPath = validationDataLocation + "exampleGATKReport.eval";
        GATKReport report = new GATKReport(reportPath);

        GATKReportTable countVariants = report.getTable("CountVariants");
        Assert.assertEquals(countVariants.getVersion(), GATKReportVersion.V0_1);
        Object countVariantsPK = countVariants.getPrimaryKey("none.eval.none.all");
        Assert.assertEquals(countVariants.get(countVariantsPK, "nProcessedLoci"), "100000");
        Assert.assertEquals(countVariants.get(countVariantsPK, "nNoCalls"), "99872");

        GATKReportTable validationReport = report.getTable("ValidationReport");
        Assert.assertEquals(validationReport.getVersion(), GATKReportVersion.V0_1);
        Object validationReportPK = countVariants.getPrimaryKey("none.eval.none.known");
        Assert.assertEquals(validationReport.get(validationReportPK, "sensitivity"), "NaN");
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

    @Test
    public void testGATKReportGatherer() {

        /*
        GATKReportTable actual1 = new GATKReportTable("TableName", "Description");
        actual1.addPrimaryKey("key");
        actual1.addColumn("colA", 0);
        actual1.addColumn("colB", 0);
        actual1.set("row1", "colA", 1);
        actual1.set("row1", "colB", 2);

        GATKReportTable actual2 = new GATKReportTable("TableName", "Description");
        actual2.addPrimaryKey("key");
        actual2.addColumn("colA", 0);
        actual2.addColumn("colB", 0);
        actual2.set("row2", "colA", 3);
        actual2.set("row2", "colB", 4);

        GATKReportTable actual3 = new GATKReportTable("TableName", "Description");
        actual3.addPrimaryKey("key");
        actual3.addColumn("colA", 0);
        actual3.addColumn("colB", 0);
        actual3.set("row3", "colA", 5);
        actual3.set("row3", "colB", 6);

        actual1.mergeRows(actual2);
        actual1.mergeRows(actual3);
        actual1.write(System.out);
        */

        GATKReportTable expected = new GATKReportTable("TableName", "Description");
        expected.addPrimaryKey("key");
        expected.addColumn("colA", 0);
        expected.addColumn("colB", 0);
        expected.set("row1", "colA", 1);
        expected.set("row1", "colB", 2);
        expected.set("row2", "colA", 3);
        expected.set("row2", "colB", 4);
        expected.set("row3", "colA", 5);
        expected.set("row3", "colB", 6);
        expected.write(System.out);

        GATKReport report1, report2, report3;
        report1 = new GATKReport();
        report1.addTable("TableName", "Description");
        report1.getTable("TableName").addPrimaryKey("key");
        report1.getTable("TableName").addColumn("colA", 0);
        report1.getTable("TableName").addColumn("colB", 0);
        report1.getTable("TableName").set("row1", "colA", 1);
        report1.getTable("TableName").set("row1", "colB", 2);

        report2 = new GATKReport();
        report2.addTable("TableName", "Description");
        report2.getTable("TableName").addPrimaryKey("key");
        report2.getTable("TableName").addColumn("colA", 0);
        report2.getTable("TableName").addColumn("colB", 0);
        report2.getTable("TableName").set("row2", "colA", 3);
        report2.getTable("TableName").set("row2", "colB", 4);

        report3 = new GATKReport();
        report3.addTable("TableName", "Description");
        report3.getTable("TableName").addPrimaryKey("key");
        report3.getTable("TableName").addColumn("colA", 0);
        report3.getTable("TableName").addColumn("colB", 0);
        report3.getTable("TableName").set("row3", "colA", 5);
        report3.getTable("TableName").set("row3", "colB", 6);

        report1.combineWith(report2);
        report1.combineWith(report3);

        report1.print(System.out);
        /*
          File a = new File("/home/roger/tbls/a.tbl");
          File b = new File("/home/roger/tbls/b.tbl");
          File c = new File("/home/roger/tbls/c.tbl");
          File out = new File("/home/roger/tbls/out.tbl");


          List<File> FileList = new ArrayList<File>();
          FileList.add(a);
          FileList.add(b);
          FileList.add(c);

          GATKReportGatherer gatherer = new GATKReportGatherer();
          gatherer.gather(FileList, out);
          System.out.print(out);
        */

        //Assert.assertEquals(1,1);

    }
}