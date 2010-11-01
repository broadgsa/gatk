package org.broadinstitute.sting.utils.report;

import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.report.tags.Param;

import org.testng.annotations.Test;


/**
 * @author aaron
 *         <p/>
 *         Class AnalysisModuleScannerUnitTest
 *         <p/>
 *         Test out the analysis scanner, which takes an analysis module and extracts out basic data
 */
public class AnalysisModuleScannerUnitTest extends BaseTest {

    @Test
    public void testBasicScan() {
        AnalysisModuleScanner scanner = new AnalysisModuleScanner(FakeAnalysis.class);

        // check we found one param, and check its description
        Assert.assertEquals(scanner.getParameters().size(), 3);
        Assert.assertTrue("basic description".equals(scanner.getParameters().values().iterator().next().description()));

        // check that the analysis name and description were set
        Assert.assertTrue("testAnalysis".equals(scanner.getAnalysis().name()));
        Assert.assertTrue("The is just a simple description".equals(scanner.getAnalysis().description()));

    }
}

// --------------------------------------------------------------------------------
// my fake analysis class
// --------------------------------------------------------------------------------
@Analysis(name = "testAnalysis", description = "The is just a simple description")
class FakeAnalysis {

    @Param(description = "basic description")
    public String text = "GRRR";

    @Param(description = "basic description")
    public String text2superlonganme = "GRRR";

    @Param(description = "basic description")
    public String text3 = "GRRR";

    @DataPoint(description = "basic description")
    public String text4 = "GRRR";

    @DataPoint(description = "basic description")
    public String text5 = "GRRR";

    @DataPoint(description = "basic description")
    public String text6 = "GRRR";

    public FakeAnalysis() {       
    }

}
