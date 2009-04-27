package org.broadinstitute.sting.utils;

import cern.jet.math.Arithmetic;

/**
 * MathUtils is a static class (no instantiation allowed!) with some useful math methods.
 *
 * @author Kiran Garimella
 */
public class MathUtils {
    /** Private constructor.  No instantiating this class! */
    private MathUtils() {}

    /**
     * Compares double values for equality (within 1e-6), or inequality.
     *
     * @param a  the first double value
     * @param b  the second double value
     * @return   -1 if a is greater than b, 0 if a is equal to be within 1e-6, 1 if b is greater than a.
     */
    public static byte compareDoubles(double a, double b) { return compareDoubles(a, b, 1e-6); }

    /**
     * Compares double values for equality (within epsilon), or inequality.
     *
     * @param a       the first double value
     * @param b       the second double value
     * @param epsilon the precision within which two double values will be considered equal
     * @return        -1 if a is greater than b, 0 if a is equal to be within epsilon, 1 if b is greater than a.
     */
    public static byte compareDoubles(double a, double b, double epsilon)
    {
        if (Math.abs(a - b) < epsilon) { return 0; }
        if (a > b) { return -1; }
        return 1;
    }

    /**
     * Compares float values for equality (within 1e-6), or inequality.
     *
     * @param a  the first float value
     * @param b  the second float value
     * @return   -1 if a is greater than b, 0 if a is equal to be within 1e-6, 1 if b is greater than a.
     */
    public static byte compareFloats(float a, float b) { return compareFloats(a, b, 1e-6f); }

    /**
     * Compares float values for equality (within epsilon), or inequality.
     *
     * @param a       the first float value
     * @param b       the second float value
     * @param epsilon the precision within which two float values will be considered equal
     * @return        -1 if a is greater than b, 0 if a is equal to be within epsilon, 1 if b is greater than a.
     */
    public static byte compareFloats(float a, float b, float epsilon)
    {
        if (Math.abs(a - b) < epsilon) { return 0; }
        if (a > b) { return -1; }
        return 1;
    }

    public static double NormalDistribution(double mean, double sd, double x)
    {
        double a = 1.0 / (sd*Math.sqrt(2.0 * Math.PI));
        double b = Math.exp(-1.0 * (Math.pow(x - mean,2.0)/(2.0 * sd * sd)));
        return a * b;
    }

    /**
     * Computes a binomial probability.  This is computed using the formula
     *
     *  B(k; n; p) = [ n! / ( k! (n - k)! ) ] (p^k)( (1-p)^k )
     *
     * where n is the number of trials, k is the number of successes, and p is the probability of success
     *
     * @param k  number of successes
     * @param n  number of Bernoulli trials
     * @param p  probability of success
     *
     * @return   the binomial probability of the specified configuration.  Computes values down to about 1e-237.
     */
    public static double binomialProbability(int k, int n, double p) {
        return Arithmetic.binomial(n, k)*Math.pow(p, k)*Math.pow(1.0 - p, n - k);
        //return (new cern.jet.random.Binomial(n, p, cern.jet.random.engine.RandomEngine.makeDefault())).pdf(k);
    }

    /**
     * Computes a multinomial probability.  This is computed using the formula
     *
     *   M(x1,x2,...,xk; n; p1,p2,...,pk) = [ n! / (x1! x2! ... xk!) ] (p1^x1)(p2^x2)(...)(pk^xk)
     *
     * where xi represents the number of times outcome i was observed, n is the number of total observations, and
     * pi represents the probability of the i-th outcome to occur.  In this implementation, the value of n is
     * inferred as the sum over i of xi;
     *
     * @param x  an int[] of counts, where each element represents the number of times a certain outcome was observed
     * @param p  a double[] of probabilities, where each element represents the probability a given outcome can occur
     * @return   the multinomial probability of the specified configuration.
     */
    public static double multinomialProbability(int[] x, double[] p) {
        int n = 0;
        for ( int obsCount : x ) { n += obsCount; }
        double nfact = Arithmetic.factorial(n);

        double obsfact = 1.0, probs = 1.0, totalprob = 0.0;
        for (int obsCountsIndex = 0; obsCountsIndex < x.length; obsCountsIndex++) {
            double ofact = Arithmetic.factorial(x[obsCountsIndex]);
            obsfact *= ofact;
            probs *= Math.pow(p[obsCountsIndex], x[obsCountsIndex]);

            totalprob += p[obsCountsIndex];
        }

        assert(MathUtils.compareDoubles(totalprob, 1.0, 0.01) == 0);

        return (nfact/obsfact)*probs;
    }


}
