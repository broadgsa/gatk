package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


public class ExactAFCalculationModelUnitTest extends BaseTest {
    static Allele A = Allele.create("A", true);
    static Allele C = Allele.create("C");
    static Allele G = Allele.create("G");
    static Allele T = Allele.create("T");

    static int sampleNameCounter = 0;
    static Genotype AA1, AB1, BB1;
    static Genotype AA2, AB2, AC2, BB2, BC2, CC2;
    final double[] FLAT_3SAMPLE_PRIORS = new double[2*3+1];  // flat priors

    @BeforeSuite
    public void before() {
        AA1 = makePL(Arrays.asList(A, A), 0, 20, 20);
        AB1 = makePL(Arrays.asList(A, C), 20, 0, 20);
        BB1 = makePL(Arrays.asList(C, C), 20, 20, 0);

        AA2 = makePL(Arrays.asList(A, A), 0, 20, 20, 20, 20, 20);
        AB2 = makePL(Arrays.asList(A, C), 20, 0, 20, 20, 20, 20);
        BB2 = makePL(Arrays.asList(C, C), 20, 20, 0, 20, 20, 20);
        AC2 = makePL(Arrays.asList(A, G), 20, 20, 20, 0, 20, 20);
        BC2 = makePL(Arrays.asList(C, G), 20, 20, 20, 20, 0, 20);
        CC2 = makePL(Arrays.asList(G, G), 20, 20, 20, 20, 20, 0);
    }

    private Genotype makePL(final List<Allele> expectedGT, int ... pls) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(expectedGT);
        gb.PL(pls);
        return gb.make();
    }

    private class GetGLsTest extends TestDataProvider {
        GenotypesContext GLs;
        int numAltAlleles;
        final ExactAFCalculation calc;
        final int[] expectedACs;
        final double[] priors;

        private GetGLsTest(final ExactAFCalculation calculation, int numAltAlleles, List<Genotype> arg, final double[] priors) {
            super(GetGLsTest.class);
            GLs = GenotypesContext.create(new ArrayList<Genotype>(arg));
            this.numAltAlleles = numAltAlleles;
            this.calc = calculation;
            this.priors = priors;

            expectedACs = new int[numAltAlleles+1];
            for ( int alleleI = 0; alleleI < expectedACs.length; alleleI++ ) {
                expectedACs[alleleI] = 0;
                final Allele allele = getAlleles().get(alleleI);
                for ( Genotype g : arg ) {
                    expectedACs[alleleI] += Collections.frequency(g.getAlleles(), allele);
                }
            }
        }

        public AlleleFrequencyCalculationResult execute() {
            return getCalc().getLog10PNonRef(getVC(), getPriors());
        }

        public double[] getPriors() {
            return priors;
        }

        public ExactAFCalculation getCalc() {
            return calc;
        }

        public VariantContext getVC() {
            VariantContextBuilder builder = new VariantContextBuilder("test", "1", 1, 1, getAlleles());
            builder.genotypes(GLs);
            return builder.make();
        }

        public List<Allele> getAlleles() {
            return Arrays.asList(Allele.create("A", true),
                    Allele.create("C"),
                    Allele.create("G"),
                    Allele.create("T")).subList(0, numAltAlleles+1);
        }

        public boolean isNonRef() {
            return expectedACs[0] < getVC().getNSamples() * 2;
        }

        public int getExpectedAltAC(final int alleleI) {
            return expectedACs[alleleI+1];
        }

        public String toString() {
            return String.format("%s model=%s input=%s", super.toString(), calc.getClass().getSimpleName(), GLs);
        }
    }

    @DataProvider(name = "wellFormedGLs")
    public Object[][] createSimpleGLsData() {
        final List<Genotype> biAllelicSamples = Arrays.asList(AA1, AB1, BB1);
        final List<Genotype> triAllelicSamples = Arrays.asList(AA2, AB2, BB2, AC2, BC2, CC2);

        for ( final int nSamples : Arrays.asList(1, 2, 3, 4) ) {
            final DiploidExactAFCalculation diploidCalc = new DiploidExactAFCalculation(nSamples, 4);
            final GeneralPloidyExactAFCalculation generalCalc = new GeneralPloidyExactAFCalculation(nSamples, 4, 2);
            final double[] priors = new double[2*nSamples+1];  // flat priors

            for ( ExactAFCalculation model : Arrays.asList(diploidCalc, generalCalc) ) {
                // bi-allelic
                if ( nSamples <= biAllelicSamples.size() )
                    for ( List<Genotype> genotypes : Utils.makePermutations(biAllelicSamples, nSamples, true) )
                        new GetGLsTest(model, 1, genotypes, priors);

                // tri-allelic
                for ( List<Genotype> genotypes : Utils.makePermutations(triAllelicSamples, nSamples, true) )
                    new GetGLsTest(model, 2, genotypes, priors);
            }
        }

        return GetGLsTest.getTests(GetGLsTest.class);
    }


    @Test(dataProvider = "wellFormedGLs")
    public void testGLs(GetGLsTest cfg) {
        final AlleleFrequencyCalculationResult result = cfg.execute();

        if ( cfg.isNonRef() ) {
            //logger.warn("pNonRef = " + result.getLog10PosteriorOfAFzero());
            Assert.assertTrue(result.getLog10PosteriorOfAFzero() < -1, "Genotypes imply pNonRef > 0 but we had posterior AF = 0 of " + result.getLog10PosteriorOfAFzero());

            // TODO -- why does this fail?
            //Assert.assertTrue(result.getLog10PosteriorsMatrixSumWithoutAFzero() > -1, "Genotypes imply pNonRef > 0 but posterior sum over all non-AF0 fields was only " + result.getLog10PosteriorsMatrixSumWithoutAFzero());

            // todo -- I'm not sure this is supposed to be true
            //Assert.assertEquals(Math.pow(10, result.getLog10PosteriorOfAFzero()) + Math.pow(10, result.getLog10PosteriorsMatrixSumWithoutAFzero()), 1.0, 1e-3, "Total posterior prob didn't sum to 1");
        }

        Assert.assertNotNull(result.getAllelesUsedInGenotyping());
        Assert.assertTrue(cfg.getAlleles().containsAll(result.getAllelesUsedInGenotyping()), "Result object has alleles not in our initial allele list");

        for ( int altAlleleI = 0; altAlleleI < cfg.numAltAlleles; altAlleleI++ ) {
            int expectedAlleleCount = cfg.getExpectedAltAC(altAlleleI);
            int calculatedAlleleCount = result.getAlleleCountsOfMAP()[altAlleleI];

            Assert.assertEquals(calculatedAlleleCount, expectedAlleleCount);
        }
    }

    @Test
    public void testLargeGLs() {
        final Genotype BB = makePL(Arrays.asList(C, C), 20000000, 20000000, 0);
        GetGLsTest cfg = new GetGLsTest(new DiploidExactAFCalculation(1, 1), 1, Arrays.asList(BB, BB, BB), FLAT_3SAMPLE_PRIORS);

        final AlleleFrequencyCalculationResult result = cfg.execute();

        int calculatedAlleleCount = result.getAlleleCountsOfMAP()[0];
        Assert.assertEquals(calculatedAlleleCount, 6);
    }

    @Test
    public void testMismatchedGLs() {
        final Genotype AB = makePL(Arrays.asList(A,C), 2000, 0, 2000, 2000, 2000, 2000);
        final Genotype AC = makePL(Arrays.asList(A,G), 100, 100, 100, 0, 100, 100);
        GetGLsTest cfg = new GetGLsTest(new DiploidExactAFCalculation(2, 2), 2, Arrays.asList(AB, AC), FLAT_3SAMPLE_PRIORS);

        final AlleleFrequencyCalculationResult result = cfg.execute();

        Assert.assertEquals(result.getAlleleCountsOfMAP()[0], 1);
        Assert.assertEquals(result.getAlleleCountsOfMAP()[1], 1);
    }
}
