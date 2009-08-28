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
        Assert.assertEquals(DiploidGenotypePriors.HUMAN_HETEROZYGOSITY,1e-3, DELTA);
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
        double[] vals = DiploidGenotypePriors.heterozygosity2DiploidProbabilities(h);
        Assert.assertEquals(vals[0], homRef, DELTA);
        Assert.assertEquals(vals[1], het, DELTA);
        Assert.assertEquals(vals[2], homVar, DELTA);
        Assert.assertEquals(DiploidGenotypePriors.heterozygosity2HomRefProbability(h), homRef, DELTA);
        Assert.assertEquals(DiploidGenotypePriors.heterozygosity2HetProbability(h), het, DELTA);
        Assert.assertEquals(DiploidGenotypePriors.heterozygosity2HomVarProbability(h), homVar, DELTA);
    }

    // 
    @Test
    public void testGenotypePriorsReferenceIndependent() {
        logger.warn("Executing testGenotypePriorsReferenceIndependent");
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
            if ( g.isHomRef(ref) ) val = DiploidGenotypePriors.heterozygosity2HomRefProbability(h);
            if ( g.isHet() )       val = DiploidGenotypePriors.heterozygosity2HetProbability(h);
            if ( g.isHomVar(ref) ) val = DiploidGenotypePriors.heterozygosity2HomVarProbability(h);

            val = log10(val);
            double e = array[g.ordinal()];
            Assert.assertEquals(String.format("%s should have p=%f but has p=%f", g, val, e), val, e, DELTA);
        }
    }

    @Test
    public void testGenotypePriorsReferencePolarized() {
        logger.warn("Executing testGenotypePriorsReferencePolarized");
        // AA, AC, AG, AT, CC, CG, CT, GG, GT, TT
        double[] array1 = {0.9985, 0.00033, 0.00033, 0.00033, 0.000166666666666667, 3.33333333333333e-06, 3.33333333333333e-06, 0.000166666666666667, 3.33333333333333e-06, 0.000166666666666667};
        logger.warn("  Array 1");
        testPolarizedGenotypePriors('A', 1e-3, 1e-5, array1);
        double[] array2 = {0.9985, 0.000333, 0.000333, 0.000333, 0.000166666666666667, 3.33333333333333e-07, 3.33333333333333e-07, 0.000166666666666667, 3.33333333333333e-07, 0.000166666666666667};
        logger.warn("  Array 2");
        testPolarizedGenotypePriors('A', 1e-3, 1e-6, array2);
        double[] array3 = {0.985, 0.00333, 0.00333, 0.00333, 0.00166666666666667, 3.33333333333333e-06, 3.33333333333333e-06, 0.00166666666666667, 3.33333333333333e-06, 0.00166666666666667};
        logger.warn("  Array 3");
        testPolarizedGenotypePriors('A', 1e-2, 1e-5, array3);
        double[] array4 = {0.99985, 3.3e-05, 3.3e-05, 3.3e-05, 1.66666666666667e-05, 3.33333333333333e-07, 3.33333333333333e-07, 1.66666666666667e-05, 3.33333333333333e-07, 1.66666666666667e-05};
        logger.warn("  Array 4");
        testPolarizedGenotypePriors('A', 1e-4, 1e-6, array4);
    }

    private void testPolarizedGenotypePriors(char ref, double h, double pTri, double[] array) {
        DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, h, pTri);
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double val = Math.pow(10, priors.getPrior(g));
            double e = array[g.ordinal()];
            Assert.assertEquals(String.format("%s should have p=%f but has p=%f", g, val, e), val, e, DELTA);
        }
    }
}