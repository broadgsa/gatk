package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class ExactAFCalculationModelUnitTest extends BaseTest {
    static Allele A = Allele.create("A", true);
    static Allele C = Allele.create("C");
    static Allele G = Allele.create("G");
    static Allele T = Allele.create("T");

    static int sampleNameCounter = 0;
    static Genotype AA1, AB1, BB1, NON_INFORMATIVE1;
    static Genotype AA2, AB2, AC2, BB2, BC2, CC2, NON_INFORMATIVE2;
    final double[] FLAT_3SAMPLE_PRIORS = new double[2*3+1];  // flat priors
    final private static boolean INCLUDE_BIALLELIC = true;
    final private static boolean INCLUDE_TRIALLELIC = true;
    final private static boolean Guillermo_FIXME = false; // TODO -- can only be enabled when GdA fixes bug

    @BeforeSuite
    public void before() {
        AA1 = makePL(Arrays.asList(A, A), 0, 20, 20);
        AB1 = makePL(Arrays.asList(A, C), 20, 0, 20);
        BB1 = makePL(Arrays.asList(C, C), 20, 20, 0);
        NON_INFORMATIVE1 = makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), 0, 0, 0);

        AA2 = makePL(Arrays.asList(A, A), 0, 20, 20, 20, 20, 20);
        AB2 = makePL(Arrays.asList(A, C), 20, 0, 20, 20, 20, 20);
        BB2 = makePL(Arrays.asList(C, C), 20, 20, 0, 20, 20, 20);
        AC2 = makePL(Arrays.asList(A, G), 20, 20, 20, 0, 20, 20);
        BC2 = makePL(Arrays.asList(C, G), 20, 20, 20, 20, 0, 20);
        CC2 = makePL(Arrays.asList(G, G), 20, 20, 20, 20, 20, 0);
        NON_INFORMATIVE2 = makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), 0, 0, 0, 0, 0, 0);
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
        final String priorName;

        private GetGLsTest(final ExactAFCalculation calculation, int numAltAlleles, List<Genotype> arg, final double[] priors, final String priorName) {
            super(GetGLsTest.class);
            GLs = GenotypesContext.create(new ArrayList<Genotype>(arg));
            this.numAltAlleles = numAltAlleles;
            this.calc = calculation;
            this.priors = priors;
            this.priorName = priorName;

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
            return String.format("%s model=%s prior=%s input=%s", super.toString(), calc.getClass().getSimpleName(),
                    priorName, GLs.size() > 5 ? String.format("%d samples", GLs.size()) : GLs);
        }
    }

    @DataProvider(name = "wellFormedGLs")
    public Object[][] createSimpleGLsData() {
        final List<Genotype> biAllelicSamples = Arrays.asList(AA1, AB1, BB1);
        final List<Genotype> triAllelicSamples = Arrays.asList(AA2, AB2, BB2, AC2, BC2, CC2);

        for ( final int nSamples : Arrays.asList(1, 2, 3, 4) ) {
            final DiploidExactAFCalculation diploidCalc = new DiploidExactAFCalculation(nSamples, 4);
            final GeneralPloidyExactAFCalculation generalCalc = new GeneralPloidyExactAFCalculation(nSamples, 4, 2);

            final int nPriorValues = 2*nSamples+1;
            final double[] flatPriors = MathUtils.normalizeFromLog10(new double[nPriorValues], true);  // flat priors
            final double[] humanPriors = new double[nPriorValues];
            UnifiedGenotyperEngine.computeAlleleFrequencyPriors(nPriorValues-1, humanPriors, 0.001);

            for ( final double[] priors : Arrays.asList(flatPriors, humanPriors) ) { // , humanPriors) ) {
                for ( ExactAFCalculation model : Arrays.asList(diploidCalc, generalCalc) ) {
                    final String priorName = priors == humanPriors ? "human" : "flat";

                    // bi-allelic
                    if ( INCLUDE_BIALLELIC && nSamples <= biAllelicSamples.size() )
                        for ( List<Genotype> genotypes : Utils.makePermutations(biAllelicSamples, nSamples, true) )
                            new GetGLsTest(model, 1, genotypes, priors, priorName);

                    // tri-allelic
                    if ( INCLUDE_TRIALLELIC && ( ! priorName.equals("human") || model != generalCalc || Guillermo_FIXME ) )
                        for ( List<Genotype> genotypes : Utils.makePermutations(triAllelicSamples, nSamples, true) )
                            new GetGLsTest(model, 2, genotypes, priors, priorName);
                }
            }
        }

        return GetGLsTest.getTests(GetGLsTest.class);
    }

    private static class NonInformativeData {
        final Genotype nonInformative;
        final List<Genotype> called;
        final int nAltAlleles;

        private NonInformativeData(List<Genotype> called, Genotype nonInformative, int nAltAlleles) {
            this.called = called;
            this.nonInformative = nonInformative;
            this.nAltAlleles = nAltAlleles;
        }
    }

    @DataProvider(name = "GLsWithNonInformative")
    public Object[][] makeGLsWithNonInformative() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<NonInformativeData> nonInformativeTests = new LinkedList<NonInformativeData>();
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB1), NON_INFORMATIVE1, 1));
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB2), NON_INFORMATIVE2, 2));
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB2, BC2), NON_INFORMATIVE2, 2));

        for ( final int nNonInformative : Arrays.asList(1, 10, 100) ) {
            for ( final NonInformativeData testData : nonInformativeTests ) {
                final List<Genotype> samples = new ArrayList<Genotype>();
                samples.addAll(testData.called);
                samples.addAll(Collections.nCopies(nNonInformative, testData.nonInformative));

                final int nSamples = samples.size();
                final DiploidExactAFCalculation diploidCalc = new DiploidExactAFCalculation(nSamples, 4);
                final GeneralPloidyExactAFCalculation generalCalc = new GeneralPloidyExactAFCalculation(nSamples, 4, 2);
                final double[] priors = new double[2*nSamples+1];  // flat priors

                for ( ExactAFCalculation model : Arrays.asList(diploidCalc, generalCalc) ) {
                    final GetGLsTest onlyInformative = new GetGLsTest(model, testData.nAltAlleles, testData.called, priors, "flat");

                    for ( int rotation = 0; rotation < nSamples; rotation++ ) {
                        Collections.rotate(samples, 1);
                        final GetGLsTest withNonInformative = new GetGLsTest(model, testData.nAltAlleles, samples, priors, "flat");
                        tests.add(new Object[]{onlyInformative, withNonInformative});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "wellFormedGLs")
    public void testGLs(GetGLsTest cfg) {
        testResultSimple(cfg);
    }

    @Test(enabled = true, dataProvider = "GLsWithNonInformative", dependsOnMethods = "testGLs")
    public void testGLsWithNonInformative(GetGLsTest onlyInformative, GetGLsTest withNonInformative) {
        final AlleleFrequencyCalculationResult expected = onlyInformative.execute();
        final AlleleFrequencyCalculationResult actual = withNonInformative.execute();

        testResultSimple(withNonInformative);

        Assert.assertEquals(actual.getLog10PosteriorOfAFzero(), expected.getLog10LikelihoodOfAFzero());
        Assert.assertEquals(actual.getLog10LikelihoodOfAFzero(), expected.getLog10LikelihoodOfAFzero());
        Assert.assertEquals(actual.getLog10PosteriorsMatrixSumWithoutAFzero(), expected.getLog10PosteriorsMatrixSumWithoutAFzero());
        Assert.assertEquals(actual.getAlleleCountsOfMAP(), expected.getAlleleCountsOfMAP());
        Assert.assertEquals(actual.getAlleleCountsOfMLE(), expected.getAlleleCountsOfMLE());
        Assert.assertEquals(actual.getLog10MAP(), expected.getLog10MAP());
        Assert.assertEquals(actual.getLog10MLE(), expected.getLog10MLE());
        Assert.assertEquals(actual.getAllelesUsedInGenotyping(), expected.getAllelesUsedInGenotyping());
    }

    private void testResultSimple(final GetGLsTest cfg) {
        final AlleleFrequencyCalculationResult result = cfg.execute();

        Assert.assertEquals(result.getNormalizedPosteriorOfAFzero() + result.getNormalizedPosteriorOfAFGTZero(), 1.0, 1e-4);

        final int minNumberOfEvaluations = cfg.getVC().getCalledChrCount();
        Assert.assertTrue(result.getnEvaluations() >= minNumberOfEvaluations,
                "Number of evaluations " + result.getnEvaluations() + " must be at least " + minNumberOfEvaluations);
        Assert.assertNotNull(result.getAllelesUsedInGenotyping());
        Assert.assertTrue(cfg.getAlleles().containsAll(result.getAllelesUsedInGenotyping()), "Result object has alleles not in our initial allele list");

        for ( int altAlleleI = 0; altAlleleI < cfg.numAltAlleles; altAlleleI++ ) {
            int expectedAlleleCount = cfg.getExpectedAltAC(altAlleleI);
            int calcAC_MLE = result.getAlleleCountsOfMLE()[altAlleleI];

            final Allele allele = cfg.getAlleles().get(altAlleleI+1);
            Assert.assertEquals(calcAC_MLE, expectedAlleleCount, "MLE AC not equal to expected AC for allele " + allele);
        }

        // TODO
        // TODO -- enable when we understand the contract between AC_MAP and pNonRef
        // TODO
//        final int AC_MAP = (int)MathUtils.sum(result.getAlleleCountsOfMAP());
//        if ( AC_MAP > 0  ) {
//            Assert.assertTrue(result.getNormalizedPosteriorOfAFzero() < 0.50, "MAP AC " + AC_MAP + " > 0 but we had posterior AF = 0 > 0.5 of " + result.getNormalizedPosteriorOfAFzero());
//        } else {
//            Assert.assertTrue(result.getNormalizedPosteriorOfAFzero() > 0.50, "MAP AC " + AC_MAP + " == 0 but we had posterior AF = 0 < 0.5 of " + result.getNormalizedPosteriorOfAFzero());
//        }
    }

    @Test(enabled = true)
    public void testLargeGLs() {
        final Genotype BB = makePL(Arrays.asList(C, C), 20000000, 20000000, 0);
        GetGLsTest cfg = new GetGLsTest(new DiploidExactAFCalculation(1, 1), 1, Arrays.asList(BB, BB, BB), FLAT_3SAMPLE_PRIORS, "flat");

        final AlleleFrequencyCalculationResult result = cfg.execute();

        int calculatedAlleleCount = result.getAlleleCountsOfMAP()[0];
        Assert.assertEquals(calculatedAlleleCount, 6);
    }

    @Test(enabled = true)
    public void testMismatchedGLs() {
        final Genotype AB = makePL(Arrays.asList(A,C), 2000, 0, 2000, 2000, 2000, 2000);
        final Genotype AC = makePL(Arrays.asList(A,G), 100, 100, 100, 0, 100, 100);
        GetGLsTest cfg = new GetGLsTest(new DiploidExactAFCalculation(2, 2), 2, Arrays.asList(AB, AC), FLAT_3SAMPLE_PRIORS, "flat");

        final AlleleFrequencyCalculationResult result = cfg.execute();

        Assert.assertEquals(result.getAlleleCountsOfMAP()[0], 1);
        Assert.assertEquals(result.getAlleleCountsOfMAP()[1], 1);
    }

    @DataProvider(name = "Models")
    public Object[][] makeModels() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new DiploidExactAFCalculation(1, 4)});
        tests.add(new Object[]{new GeneralPloidyExactAFCalculation(1, 4, 2)});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "Models")
    public void testBiallelicPriors(final ExactAFCalculation model) {
        final int REF_PL = 10;
        final Genotype AB = makePL(Arrays.asList(A,C), REF_PL, 0, 10000);

        for ( int log10NonRefPrior = 1; log10NonRefPrior < 100*REF_PL; log10NonRefPrior += 1 ) {
            final double refPrior = 1 - QualityUtils.qualToErrorProb(log10NonRefPrior);
            final double[] priors = MathUtils.toLog10(new double[]{refPrior, (1-refPrior) / 2, (1-refPrior) / 2});
            GetGLsTest cfg = new GetGLsTest(model, 1, Arrays.asList(AB), priors, "pNonRef" + log10NonRefPrior);
            final AlleleFrequencyCalculationResult result = cfg.execute();
            final int actualAC = result.getAlleleCountsOfMAP()[0];

            final double pRefWithPrior = AB.getLikelihoods().getAsVector()[0] + priors[0];
            final double pHetWithPrior = AB.getLikelihoods().getAsVector()[1] + priors[1];
            final boolean expectNonRef = pRefWithPrior <= pHetWithPrior;

            if ( expectNonRef )
                Assert.assertTrue(result.getNormalizedPosteriorOfAFGTZero() > 0.5);
            else
                Assert.assertTrue(result.getNormalizedPosteriorOfAFGTZero() < 0.5);

            final int expectedAC = expectNonRef ? 1 : 0;
            Assert.assertEquals(actualAC, expectedAC,
                    "actual AC with priors " + log10NonRefPrior + " not expected "
                            + expectedAC + " priors " + Utils.join(",", priors));
        }
    }

    @Test(enabled = false, dataProvider = "Models")
    public void testTriallelicPriors(final ExactAFCalculation model) {
        // TODO
        // TODO
        // TODO THIS SEEMS TO ID A BUG IN THE EXACT MODEL FOR MULTI-ALLELICS, AS THE
        // TODO SECOND ALLELE ISN'T HAVING A SQUARED PRIOR.  TALK TO ERIC AND CONFIRM
        // TODO
        // TODO
        final int REF_PL_AB = 10, REF_PL_AC = 20; // first AC goes, then AB
        final Genotype AB = makePL(Arrays.asList(A,C), REF_PL_AB, 0, 10000, 10000, 10000);
        final Genotype AC = makePL(Arrays.asList(A, G), REF_PL_AC, 10000, 10000, 0, 10000, 10000);

        for ( int log10NonRefPrior = 1; log10NonRefPrior < 100*REF_PL_AC; log10NonRefPrior += 1 ) {
            final double refPrior = 1 - QualityUtils.qualToErrorProb(log10NonRefPrior);
            final double nonRefPrior = (1-refPrior) / 2;
            final double[] priors = MathUtils.toLog10(new double[]{refPrior, nonRefPrior, nonRefPrior, nonRefPrior, nonRefPrior, nonRefPrior});
            GetGLsTest cfg = new GetGLsTest(model, 2, Arrays.asList(AB, AC), priors, "pNonRef" + log10NonRefPrior);
            final AlleleFrequencyCalculationResult result = cfg.execute();
            final int actualAC_AB = result.getAlleleCountsOfMAP()[0];

            final double pRefABWithPrior = AB.getLikelihoods().getAsVector()[0] + priors[0];
            final double pHetABWithPrior = AB.getLikelihoods().getAsVector()[1] + priors[1];
            final int expectedAC_AB = pRefABWithPrior <= pHetABWithPrior ? 1 : 0;
            Assert.assertEquals(actualAC_AB, expectedAC_AB,
                    "actual AC with priors " + log10NonRefPrior + " not expected "
                            + expectedAC_AB + " priors " + Utils.join(",", priors));

            final double nonRefPriorSecondAllele = Math.pow(nonRefPrior, 2);
            final double refPriorSecondAllele = 1 - nonRefPriorSecondAllele;
            final int actualAC_AC = result.getAlleleCountsOfMAP()[1];
            final double pRefACWithPrior = AB.getLikelihoods().getAsVector()[0] + Math.log10(refPriorSecondAllele);
            final double pHetACWithPrior = AC.getLikelihoods().getAsVector()[3] + Math.log10(nonRefPriorSecondAllele);
            final int expectedAC_AC = pRefACWithPrior <= pHetACWithPrior ? 1 : 0;
            Assert.assertEquals(actualAC_AC, expectedAC_AC,
                    "actual AC with priors " + log10NonRefPrior + " not expected "
                            + expectedAC_AC + " priors " + Utils.join(",", priors));
        }
    }
}