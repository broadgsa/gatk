package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalData;
import org.junit.Before;
import org.junit.Test;
import org.junit.Assert;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

public class GenotypeLikelihoodsTest extends BaseTest {
    private final static double DELTA = 1e-8;

    @Test
    public void testBasic() {
        logger.warn("Executing testIsBetween");
        Assert.assertTrue(GenotypeLikelihoods.oneMinusData.length >= Byte.MAX_VALUE);
        Assert.assertTrue(GenotypeLikelihoods.oneHalfMinusDataArachne.length >= Byte.MAX_VALUE);
        Assert.assertTrue(GenotypeLikelihoods.oneHalfMinusData3Base.length >= Byte.MAX_VALUE);
        Assert.assertEquals(GenotypeLikelihoods.log10Of1_3,-0.4771212547196624, DELTA);
        Assert.assertEquals(GenotypeLikelihoods.HUMAN_HETEROZYGOSITY,1e-3, DELTA);

        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            double e = pow(10, (qual / -10.0));
            Assert.assertEquals(GenotypeLikelihoods.oneMinusData[qual], log10(1.0 - e), DELTA);
            Assert.assertEquals(GenotypeLikelihoods.oneHalfMinusDataArachne[qual], log10(0.5 - e / 2.0), DELTA);
            Assert.assertEquals(GenotypeLikelihoods.oneHalfMinusData3Base[qual], log10(0.5 - e / 2.0 + e / 6.0), DELTA);
        }
    }


    // f <- function(h) { print(paste(1-3.0 * h / 2, h, h/2, sep=', '));}
    @Test
    public void testPriorsFromHet() {
        logger.warn("Executing testPriorsFromHet");
        testPriorsFromHet(0.0, 1, 0, 0);
        testPriorsFromHet(1e-1, 0.85, 0.1, 0.05);
        testPriorsFromHet(1e-2, 0.985, 0.01, 0.005);
        testPriorsFromHet(1e-3, 0.9985, 0.001, 5e-04);
        testPriorsFromHet(1e-4, 0.99985, 1e-04, 5e-05);
        testPriorsFromHet(1e-5, 0.999985, 1e-05, 5e-06);
        testPriorsFromHet(0.5, 0.25, 0.5, 0.25);
    }

    @Test (expected = RuntimeException.class)
    public void testPriorsFromHetFail1() {
        logger.warn("Executing testPriorsFromHetFail1");
        testPriorsFromHet(1.0, 0, 0, 0);
    }

    @Test (expected = RuntimeException.class)
    public void testPriorsFromHetFail2() {
        logger.warn("Executing testPriorsFromHetFail2");
        testPriorsFromHet(-1.0, 0, 0, 0);
    }

    private void testPriorsFromHet(double h, double homRef, double het, double homVar) {
        double[] vals = GenotypeLikelihoods.heterozygosity2DiploidProbabilities(h);
        Assert.assertEquals(vals[0], homRef, DELTA);
        Assert.assertEquals(vals[1], het, DELTA);
        Assert.assertEquals(vals[2], homVar, DELTA);
        Assert.assertEquals(GenotypeLikelihoods.heterozygosity2HomRefProbability(h), homRef, DELTA);
        Assert.assertEquals(GenotypeLikelihoods.heterozygosity2HetProbability(h), het, DELTA);
        Assert.assertEquals(GenotypeLikelihoods.heterozygosity2HomVarProbability(h), homVar, DELTA);
    }

    // 
    @Test
    public void testGenotypePriors1() {
        logger.warn("Executing testGenotypePriors1");
        // AA, AC, AG, AT, CC, CG, CT, GG, GT, TT
        double[] array1 = {-0.0705810742857073, -1, -1, -1, -1.301029995663981, -1, -1, -1.301029995663981, -1, -1.301029995663981};
        testGenotypePriors('A', 1e-1, array1);
        double[] array2 = {-1.301029995663981, -1, -1, -1, -0.0705810742857073, -1, -1, -1.301029995663981, -1, -1.301029995663981};
        testGenotypePriors('C', 1e-1, array2);
        double[] array3 = {-1.301029995663981, -1, -1, -1, -1.301029995663981, -1, -1, -0.0705810742857073, -1, -1.301029995663981};
        testGenotypePriors('G', 1e-1, array3);
        double[] array4 = {-1.301029995663981, -1, -1, -1, -1.301029995663981, -1, -1, -1.301029995663981, -1, -0.0705810742857073};
        testGenotypePriors('T', 1e-1, array4);
    }

    private void testGenotypePriors(char ref, double h, double[] array) {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double val = 0.0;
            if ( g.isHomRef(ref) ) val = GenotypeLikelihoods.heterozygosity2HomRefProbability(h);
            if ( g.isHet() )       val = GenotypeLikelihoods.heterozygosity2HetProbability(h);
            if ( g.isHomVar(ref) ) val = GenotypeLikelihoods.heterozygosity2HomVarProbability(h);

            val = log10(val);
            double e = array[g.ordinal()];
            Assert.assertEquals(String.format("%s should have p=%f but has p=%f", g, val, e), val, e, DELTA);
        }
    }


    /**
     * Takes reference base, and three priors for hom-ref, het, hom-var, and fills in the priors vector
     * appropriately.
     *
     * @param ref
     * @param priorHomRef
     * @param priorHet
     * @param priorHomVar
     */
    public static double[] getGenotypePriors(char ref, double priorHomRef, double priorHet, double priorHomVar) {
        if ((priorHomRef + priorHet + priorHomVar) != 1) {
            throw new RuntimeException(String.format("Prior probabilities don't sum to one => %f, %f, %f", priorHomRef, priorHet, priorHomVar));
        }

        double[] priors = new double[DiploidGenotype.values().length];

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int nRefBases = Utils.countOccurrences(ref, g.toString());
            double log10POfG = 0.0;

            switch ( nRefBases ) {
                case 2: // hom-ref
                    log10POfG = Math.log10(priorHomRef);
                    break;
                case 0: // hom-nonref
                    //log10POfG = Math.log10(priorHomVar / 3);
                    log10POfG = Math.log10(priorHomVar);
                    break;
                default:
                    //log10POfG = Math.log10(priorHet / 6);
                    log10POfG = Math.log10(priorHet);
                    break;
            }

            priors[g.ordinal()] = log10POfG;
        }

        return priors;
    }
}