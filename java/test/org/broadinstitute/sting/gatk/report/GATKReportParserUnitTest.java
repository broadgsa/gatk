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

import java.io.File;

public class GATKReportParserUnitTest extends BaseTest {
    @Test
    public void testParse() throws Exception {
        GATKReportParser parser = new GATKReportParser();
        parser.parse(new File(validationDataLocation + "exampleGATKReport.eval"));

        Assert.assertEquals(parser.getValue("CountVariants", "none.eval.none.all", "nProcessedLoci"), "100000");
        Assert.assertEquals(parser.getValue("CountVariants", "none.eval.none.all", "nNoCalls"), "99872");

        Assert.assertEquals(parser.getValue("SimpleMetricsByAC.metrics", "none.eval.none.novel.ac2", "AC"), "2");
        Assert.assertNull(parser.getValue("SimpleMetricsByAC.metrics", "none.eval.none.novel.ac2.bad", "AC"));
        Assert.assertNull(parser.getValue("SimpleMetricsByAC.metrics", "none.eval.none.novel.ac2", "AC.bad"));
        Assert.assertNull(parser.getValue("SimpleMetricsByAC.metrics.bad", "none.eval.none.novel.ac2", "AC"));

        Assert.assertEquals(parser.getValue("ValidationReport", "none.eval.none.known", "sensitivity"), "NaN");
    }
}
