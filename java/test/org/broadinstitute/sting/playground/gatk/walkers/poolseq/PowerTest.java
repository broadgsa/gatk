package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.junit.Test;
import org.junit.Assert;

import java.math.BigDecimal;


/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 14, 2009
 * Time: 10:19:53 AM
 * To change this template use File | Settings | File Templates.
 */
public class PowerTest extends BaseTest {
    // test for various items involved in PowerBelowFrequency
    @Test
    public void testBinomialProbabilityLog() {
        logger.warn("Testing binomial probability in log space");
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbabilityLog(5,7,0.3),0.02500,1e-4)==0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbabilityLog(8000,15000,0.6),2.56e-62,1e-64) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbabilityLog(12000,20000,0.5),7.34e-178,1e-180) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbabilityLog(600,800,0.75),0.0326,0.0001) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbabilityLog(50,500,0.2),6.404e-10,1e-13)==0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbabilityLog(5277,8279,0.60344),1.6023145e-11,1e-18) == 0);
    }

    @Test
    public void testBinomialCDF() {
        logger.warn("Testing binomial CDF");
        double test1 = MathUtils.cumBinomialProbLog(0,5277,8279,0.64344);
        System.out.println("Binomial CDF -- Test 1: "+test1);
        Assert.assertTrue(MathUtils.compareDoubles(test1,0.127,5e-3) == 0);
        double test2 = MathUtils.cumBinomialProbLog(0,1472,9886,0.172);
        System.out.println("Binomial CDF -- Test 2: "+test2);
        Assert.assertTrue(MathUtils.compareDoubles(test2,3.1e-10,5e-11)==0);
        double test3 = MathUtils.cumBinomialProbLog(0,1472,9886,0.142);
        System.out.println("Binomial CDF -- Test 3: "+test3);
        Assert.assertTrue(MathUtils.compareDoubles(test3,0.975,5e-3)==0);
        double test4 = MathUtils.cumBinomialProbLog(0,6781,21297,0.321742);
        System.out.println("Binomial CDF -- Test 4:"+test4);
        Assert.assertTrue(MathUtils.compareDoubles(test4,0.150,5e-3)==0);
        double test5 = MathUtils.cumBinomialProbLog(0,10609,21297,0.5001);
        System.out.println("Binomial CDF == Test 5:"+test5);
        Assert.assertTrue(MathUtils.compareDoubles(test5,0.286,5e-3)==0);
    }



/*    @Test
    public void testPowerCalculation() {
        PowerBelowFrequencyWalker w = new PowerBelowFrequencyWalker();
        w.initialize();
        logger.warn("Testing power on pool of 100 people");
        w.setPoolSize(100);
        // TEST 1
        double error = 0.000123;
        byte q = QualityUtils.probToQual(1-error);
        double pow = w.theoreticalPower(1000,q,1,5);
        System.out.println("Qual for "+error+"="+q+" and power="+pow);
        Assert.assertTrue(MathUtils.compareDoubles(pow,0.7356,1e-4) == 0);
        // TEST 2
        error = 0.00512;
        q = QualityUtils.probToQual(1-error);
        double err2 = QualityUtils.qualToErrorProb(q);
        pow = w.theoreticalPower(40,q,12,3);
        System.out.println("Qual for "+error+"="+q+" back to error="+err2+" and power="+pow);
        Assert.assertTrue(MathUtils.compareDoubles(pow,0.2173,1e-4) == 0);
        
    }*/
}
