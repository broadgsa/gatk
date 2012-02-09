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
    @Test(enabled = false)
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
        return new Object[][] {
                new Object[] {null, true},
                new Object[] {"null", true},
                new Object[] {"NA", true},
                new Object[] {"0", true},
                new Object[] {"0.0", true},
                new Object[] {"-0", true},
                new Object[] {"-0.0", true},
                new Object[] {String.valueOf(Long.MAX_VALUE), true},
                new Object[] {String.valueOf(Long.MIN_VALUE), true},
                new Object[] {String.valueOf(Float.MIN_NORMAL), true},
                new Object[] {String.valueOf(Double.MAX_VALUE), true},
                new Object[] {String.valueOf(Double.MIN_VALUE), true},
                new Object[] {String.valueOf(Double.POSITIVE_INFINITY), true},
                new Object[] {String.valueOf(Double.NEGATIVE_INFINITY), true},
                new Object[] {String.valueOf(Double.NaN), true},
                new Object[] {"hello", false}
        };
    }

    @Test(dataProvider = "rightAlignValues")
    public void testIsRightAlign(String value, boolean expected) {
        Assert.assertEquals(GATKReportColumn.isRightAlign(value), expected, "right align of '" + value + "'");
    }
}
