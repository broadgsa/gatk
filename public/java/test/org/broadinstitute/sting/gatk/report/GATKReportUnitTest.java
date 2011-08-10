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

        GATKReportTable simpleMetricsByAC = report.getTable("SimpleMetricsByAC.metrics");
        Assert.assertEquals(simpleMetricsByAC.getVersion(), GATKReportVersion.V0_1);
        Object simpleMetricsByACPK = simpleMetricsByAC.getPrimaryKey("none.eval.none.novel.ac2");
        Assert.assertEquals(simpleMetricsByAC.get(simpleMetricsByACPK, "AC"), "2");

        Assert.assertFalse(simpleMetricsByAC.containsPrimaryKey("none.eval.none.novel.ac2.bad"));
    }
}
