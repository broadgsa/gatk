package org.broadinstitute.sting.playground.utils.report;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.tags.Param;
import org.broadinstitute.sting.playground.utils.report.tags.Table;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class AnalysisModuleScannerTest
 *         <p/>
 *         Test out the analysis scanner, which takes an analysis module and extracts out basic data
 */
public class AnalysisModuleScannerTest extends BaseTest {

    @Test
    public void testBasicScan() {
        AnalysisModuleScanner scanner = new AnalysisModuleScanner(FakeAnalysis.class);

        // check we found one param, and check its description
        Assert.assertEquals(3, scanner.getParameters().size());
        Assert.assertTrue("basic description".equals(scanner.getParameters().values().iterator().next().description()));

        // check that we've found a table, and check its description
        Assert.assertEquals(1, scanner.getTables().size());
        Assert.assertTrue("Generate a table from this data".equals(scanner.getTables().values().iterator().next().description()));

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

    @Table(name="FakeTable", description = "Generate a table from this data",columns=3)
    public List<Integer> values = new ArrayList<Integer>();

    public FakeAnalysis() {
        values.add(1);
        values.add(2);
        values.add(3);
        values.add(4);
        values.add(5);
        values.add(6);

    }

}
