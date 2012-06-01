package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;


public class ExactAFCalculationModelUnitTest extends BaseTest {

    static double[] AA1, AB1, BB1;
    static double[] AA2, AB2, AC2, BB2, BC2, CC2;
    static final int numSamples = 3;
    static double[] priors = new double[2*numSamples+1];  // flat priors

    @BeforeSuite
    public void before() {
        AA1 = new double[]{0.0, -20.0, -20.0};
        AB1 = new double[]{-20.0, 0.0, -20.0};
        BB1 = new double[]{-20.0, -20.0, 0.0};
        AA2 = new double[]{0.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        AB2 = new double[]{-20.0, 0.0, -20.0, -20.0, -20.0, -20.0};
        AC2 = new double[]{-20.0, -20.0, -20.0, 0.0, -20.0, -20.0};
        BB2 = new double[]{-20.0, -20.0, 0.0, -20.0, -20.0, -20.0};
        BC2 = new double[]{-20.0, -20.0, -20.0, -20.0, 0.0, -20.0};
        CC2 = new double[]{-20.0, -20.0, -20.0, -20.0, -20.0, 0.0};
    }

    private class GetGLsTest extends TestDataProvider {
        GenotypesContext GLs;
        int numAltAlleles;
        String name;

        private GetGLsTest(String name, int numAltAlleles, Genotype... arg) {
            super(GetGLsTest.class, name);
            GLs = GenotypesContext.create(arg);
            this.name = name;
            this.numAltAlleles = numAltAlleles;
        }

        public String toString() {
            return String.format("%s input=%s", super.toString(), GLs);
        }
    }

    private static Genotype createGenotype(String name, double[] gls) {
        return new GenotypeBuilder(name, Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).PL(gls).make();
    }

    @DataProvider(name = "getGLs")
    public Object[][] createGLsData() {

        // bi-allelic case
        new GetGLsTest("B0", 1, createGenotype("AA1", AA1), createGenotype("AA2", AA1), createGenotype("AA3", AA1));
        new GetGLsTest("B1", 1, createGenotype("AA1", AA1), createGenotype("AA2", AA1), createGenotype("AB", AB1));
        new GetGLsTest("B2", 1, createGenotype("AA1", AA1), createGenotype("BB", BB1), createGenotype("AA2", AA1));
        new GetGLsTest("B3a", 1, createGenotype("AB", AB1), createGenotype("AA", AA1), createGenotype("BB", BB1));
        new GetGLsTest("B3b", 1, createGenotype("AB1", AB1), createGenotype("AB2", AB1), createGenotype("AB3", AB1));
        new GetGLsTest("B4", 1, createGenotype("BB1", BB1), createGenotype("BB2", BB1), createGenotype("AA", AA1));
        new GetGLsTest("B5", 1, createGenotype("BB1", BB1), createGenotype("AB", AB1), createGenotype("BB2", BB1));
        new GetGLsTest("B6", 1, createGenotype("BB1", BB1), createGenotype("BB2", BB1), createGenotype("BB3", BB1));

        // tri-allelic case
        new GetGLsTest("B1C0", 2, createGenotype("AA1", AA2), createGenotype("AA2", AA2), createGenotype("AB", AB2));
        new GetGLsTest("B0C1", 2, createGenotype("AA1", AA2), createGenotype("AA2", AA2), createGenotype("AC", AC2));
        new GetGLsTest("B1C1a", 2, createGenotype("AA", AA2), createGenotype("AB", AB2), createGenotype("AC", AC2));
        new GetGLsTest("B1C1b", 2, createGenotype("AA1", AA2), createGenotype("AA2", AA2), createGenotype("BC", BC2));
        new GetGLsTest("B2C1", 2, createGenotype("AB1", AB2), createGenotype("AB2", AB2), createGenotype("AC", AC2));
        new GetGLsTest("B3C2a", 2, createGenotype("AB", AB2), createGenotype("BC1", BC2), createGenotype("BC2", BC2));
        new GetGLsTest("B3C2b", 2, createGenotype("AB", AB2), createGenotype("BB", BB2), createGenotype("CC", CC2));

        return GetGLsTest.getTests(GetGLsTest.class);
    }


    @Test(dataProvider = "getGLs")
    public void testGLs(GetGLsTest cfg) {

        final AlleleFrequencyCalculationResult result = new AlleleFrequencyCalculationResult(2);

        ExactAFCalculationModel.linearExactMultiAllelic(cfg.GLs, cfg.numAltAlleles, priors, result);

        int nameIndex = 1;
        for ( int allele = 0; allele < cfg.numAltAlleles; allele++, nameIndex+=2 ) {
            int expectedAlleleCount = Integer.valueOf(cfg.name.substring(nameIndex, nameIndex+1));
            int calculatedAlleleCount = result.getAlleleCountsOfMAP()[allele];

            Assert.assertEquals(calculatedAlleleCount, expectedAlleleCount);
        }
    }

    @Test
    public void testLargeGLs() {

        final double[] BB = new double[]{-20000000.0, -20000000.0, 0.0};
        GetGLsTest cfg = new GetGLsTest("B6", 1, createGenotype("1", BB), createGenotype("2", BB), createGenotype("3", BB));

        final AlleleFrequencyCalculationResult result = new AlleleFrequencyCalculationResult(2);

        ExactAFCalculationModel.linearExactMultiAllelic(cfg.GLs, cfg.numAltAlleles, priors, result);

        int calculatedAlleleCount = result.getAlleleCountsOfMAP()[0];
        Assert.assertEquals(calculatedAlleleCount, 6);
    }
}
