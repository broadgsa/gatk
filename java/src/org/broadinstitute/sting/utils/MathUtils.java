/*
 * Copyright (c) 2010 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import cern.jet.math.Arithmetic;

import java.math.BigDecimal;
import java.util.*;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

/**
 * MathUtils is a static class (no instantiation allowed!) with some useful math methods.
 *
 * @author Kiran Garimella
 */
public class MathUtils {
    /** Public constants - used for the Lanczos approximation to the factorial function
     *  (for the calculation of the binomial/multinomial probability in logspace)
     * @param LANC_SEQ[] - an array holding the constants which correspond to the product
     * of Chebyshev Polynomial coefficients, and points on the Gamma function (for interpolation)
     * [see A Precision Approximation of the Gamma Function J. SIAM Numer. Anal. Ser. B, Vol. 1 1964. pp. 86-96]
     * @param LANC_G - a value for the Lanczos approximation to the gamma function that works to
     * high precision 
     */


    /** Private constructor.  No instantiating this class! */
    private MathUtils() {}

    public static double sum(Collection<Number> numbers) {
        return sum(numbers,false);
    }

    public static double sum( Collection<Number> numbers, boolean ignoreNan ) {
        double sum = 0;
        for ( Number n : numbers ) {
            if ( ! ignoreNan || ! Double.isNaN(n.doubleValue())) {
                sum += n.doubleValue();
            }
        }

        return sum;
    }

    public static int nonNanSize(Collection<Number> numbers) {
        int size = 0;
        for ( Number n : numbers) {
            size += Double.isNaN(n.doubleValue()) ? 0 : 1;
        }

        return size;
    }

    public static double average( Collection<Number> numbers, boolean ignoreNan) {
        if ( ignoreNan ) {
            return sum(numbers,true)/nonNanSize(numbers);
        } else {
            return sum(numbers,false)/nonNanSize(numbers);
        }
    }

    public static double variance( Collection<Number> numbers, Number mean, boolean ignoreNan ) {
        double mn = mean.doubleValue();
        double var = 0;
        for ( Number n : numbers ) { var += ( ! ignoreNan || ! Double.isNaN(n.doubleValue())) ? (n.doubleValue()-mn)*(n.doubleValue()-mn) : 0; }
        if ( ignoreNan ) { return var/(nonNanSize(numbers)-1); }
        return var/(numbers.size()-1);
    }

    public static double variance(Collection<Number> numbers, Number mean) {
        return variance(numbers,mean,false);
    }

    public static double variance(Collection<Number> numbers, boolean ignoreNan) {
        return variance(numbers,average(numbers,ignoreNan),ignoreNan);
    }

    public static double variance(Collection<Number> numbers) {
        return variance(numbers,average(numbers,false),false);
    }

    public static double sum(double[] values) {
        double s = 0.0;
        for ( double v : values) s += v;
        return s;
    }


    /**
     * Converts a real space array of probabilities into a log10 array
     * @param prRealSpace
     * @return
     */
    public static double[] toLog10(double[] prRealSpace) {
        double[] log10s = new double[prRealSpace.length];
        for ( int i = 0; i < prRealSpace.length; i++ )
            log10s[i] = Math.log10(prRealSpace[i]);
        return log10s;
    }

    public static double log10sumLog10(double[] log10p, int start) {
        double sum = 0.0;

        double maxValue = Utils.findMaxEntry(log10p);
        for ( int i = start; i < log10p.length; i++ ) {
            sum += Math.pow(10.0, log10p[i] - maxValue);
        }

        return Math.log10(sum) + maxValue;
    }

    public static double sum(List<Double> values) {
        double s = 0.0;
        for ( double v : values) s += v;
        return s;
    }

    public static int sum(List<Integer> values) {
        int s = 0;
        for ( int v : values) s += v;
        return s;
    }

    public static double sumLog10(double[] log10values) {
        return Math.pow(10.0, log10sumLog10(log10values));
//        double s = 0.0;
//        for ( double v : log10values) s += Math.pow(10.0, v);
//        return s;
    }

    public static double log10sumLog10(double[] log10values) {
        return log10sumLog10(log10values, 0);
    }

    public static boolean wellFormedDouble(double val) {
        return ! Double.isInfinite(val) && ! Double.isNaN(val);
    }

    public static boolean isBounded(double val, double lower, double upper) {
        return val >= lower && val <= upper;
    }

    public static boolean isPositive(double val) {
        return ! isNegativeOrZero(val);
    }

    public static boolean isPositiveOrZero(double val) {
        return isBounded(val, 0.0, Double.POSITIVE_INFINITY);
    }

    public static boolean isNegativeOrZero(double val) {
        return isBounded(val, Double.NEGATIVE_INFINITY, 0.0);
    }

    public static boolean isNegative(double val) {
        return ! isPositiveOrZero(val);
    }

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
     * Performs the calculation for a binomial probability in logspace, preventing the blowups of the binomial coefficient
     * consistent with moderately large n and k.
     *
     * @param k - number of successes/observations
     * @param n - number of Bernoulli trials
     * @param p - probability of success/observations
     *
     * @return    the probability mass associated with a binomial mass function for the above configuration
     */
    public static double binomialProbabilityLog(int k, int n, double p) {

        double log_coef = 0.0;
        int min;
        int max;

        if(k < n - k) {
            min = k;
            max = n-k;
        } else {
            min = n - k;
            max = k;
        }

        for(int i=2; i <= min; i++) {
            log_coef -= Math.log((double)i);
        }

        for(int i = max+1; i <= n; i++) {
            log_coef += Math.log((double)i);
        }

        return Math.exp(log_coef + ((double)k)*Math.log(p) + ((double)(n-k))*Math.log(1-p));

        // in the future, we may want to precompile a table of exact values, where the entry (a,b) indicates
        // the sum of log(i) from i = (1+(50*a)) to i = 50+(50*b) to increase performance on binomial coefficients
        // which may require many sums.
    }

    /**
     * Performs the cumulative sum of binomial probabilities, where the probability calculation is done in log space.
     * @param start - start of the cumulant sum (over hits)
     * @param end - end of the cumulant sum (over hits)
     * @param total - number of attempts for the number of hits
     * @param probHit - probability of a successful hit
     * @return - returns the cumulative probability
     */
    public static double cumBinomialProbLog(int start, int end, int total, double probHit) {
        double cumProb = 0.0;
        double prevProb;
        BigDecimal probCache = BigDecimal.ZERO;

        for(int hits = start; hits < end; hits++) {
            prevProb = cumProb;
            double probability = binomialProbabilityLog(hits,total,probHit);
            cumProb += probability;
            if ( probability > 0 && cumProb - prevProb < probability/2 ) { // loss of precision
                probCache = probCache.add(new BigDecimal(prevProb));
                cumProb = 0.0;
                hits--; // repeat loop
                // prevProb changes at start of loop
            }
        }

        return probCache.add(new BigDecimal(cumProb)).doubleValue();
    }
  
    /**
     * Computes a multinomial.  This is computed using the formula
     *
     *   M(x1,x2,...,xk; n) = [ n! / (x1! x2! ... xk!) ]
     *
     * where xi represents the number of times outcome i was observed, n is the number of total observations.
     * In this implementation, the value of n is inferred as the sum over i of xi.
     *
     * @param x  an int[] of counts, where each element represents the number of times a certain outcome was observed
     * @return   the multinomial of the specified configuration.
     */
    public static double multinomial(int[] x) {
        // In order to avoid overflow in computing large factorials in the multinomial
        // coefficient, we split the calculation up into the product of a bunch of
        // binomial coefficients.

        double multinomialCoefficient = 1.0;

        for (int i = 0; i < x.length; i++) {
            int n = 0;
            for (int j = 0; j <= i; j++) { n += x[j]; }

            double multinomialTerm = Arithmetic.binomial(n, x[i]);
            multinomialCoefficient *= multinomialTerm;
        }

        return multinomialCoefficient;
    }

    /**
     * Computes a multinomial probability.  This is computed using the formula
     *
     *   M(x1,x2,...,xk; n; p1,p2,...,pk) = [ n! / (x1! x2! ... xk!) ] (p1^x1)(p2^x2)(...)(pk^xk)
     *
     * where xi represents the number of times outcome i was observed, n is the number of total observations, and
     * pi represents the probability of the i-th outcome to occur.  In this implementation, the value of n is
     * inferred as the sum over i of xi.
     *
     * @param x  an int[] of counts, where each element represents the number of times a certain outcome was observed
     * @param p  a double[] of probabilities, where each element represents the probability a given outcome can occur
     * @return   the multinomial probability of the specified configuration.
     */
    public static double multinomialProbability(int[] x, double[] p) {
        // In order to avoid overflow in computing large factorials in the multinomial
        // coefficient, we split the calculation up into the product of a bunch of
        // binomial coefficients.
        double multinomialCoefficient = multinomial(x);

        double probs = 1.0, totalprob = 0.0;
        for (int obsCountsIndex = 0; obsCountsIndex < x.length; obsCountsIndex++) {
            probs *= Math.pow(p[obsCountsIndex], x[obsCountsIndex]);
            totalprob += p[obsCountsIndex];
        }

        assert(MathUtils.compareDoubles(totalprob, 1.0, 0.01) == 0);

        return multinomialCoefficient*probs;
    }

    /**
     * calculate the Root Mean Square of an array of integers
     * @param x  an int[] of numbers
     * @return   the RMS of the specified numbers.
    */
    public static double rms(int[] x) {
        if ( x.length == 0 )
            return 0.0;

        double rms = 0.0;
        for (int i : x)
            rms += i * i;
        rms /= x.length;
        return Math.sqrt(rms);
    }

    /**
     * calculate the Root Mean Square of an array of doubles
     * @param x  a double[] of numbers
     * @return   the RMS of the specified numbers.
    */
    public static double rms(Double[] x) {
        if ( x.length == 0 )
            return 0.0;

        double rms = 0.0;
        for (Double i : x)
            rms += i * i;
        rms /= x.length;
        return Math.sqrt(rms);
    }

    public static double distanceSquared( final double[] x, final double[] y ) {
        double dist = 0.0;
        for(int iii = 0; iii < x.length; iii++) {
            dist += (x[iii] - y[iii]) * (x[iii] - y[iii]);
        }
        return dist;
    }

    /**
     * normalizes the log10-based array.  ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0 (<= 1 IN REAL-SPACE).
     *
     * @param array  the array to be normalized
     * @param takeLog10OfOutput if true, the output will be transformed back into log10 units
     *
     * @return a newly allocated array corresponding the normalized values in array, maybe log10 transformed
    */
    public static double[] normalizeFromLog10(double[] array, boolean takeLog10OfOutput) {
        double[] normalized = new double[array.length];

        // for precision purposes, we need to add (or really subtract, since they're
        // all negative) the largest value; also, we need to convert to normal-space.
        double maxValue = Utils.findMaxEntry(array);
        for (int i = 0; i < array.length; i++)
            normalized[i] = Math.pow(10, array[i] - maxValue);

        // normalize
        double sum = 0.0;
        for (int i = 0; i < array.length; i++)
            sum += normalized[i];
        for (int i = 0; i < array.length; i++) {
            double x = normalized[i] / sum;
            if ( takeLog10OfOutput ) x = Math.log10(x);
            normalized[i] = x;
        }

        return normalized;
    }

    public static double[] normalizeFromLog10(List<Double> array, boolean takeLog10OfOutput) {
        double[] normalized = new double[array.size()];

        // for precision purposes, we need to add (or really subtract, since they're
        // all negative) the largest value; also, we need to convert to normal-space.
        double maxValue = MathUtils.arrayMax( array );
        for (int i = 0; i < array.size(); i++)
            normalized[i] = Math.pow(10, array.get(i) - maxValue);

        // normalize
        double sum = 0.0;
        for (int i = 0; i < array.size(); i++)
            sum += normalized[i];
        for (int i = 0; i < array.size(); i++) {
            double x = normalized[i] / sum;
            if ( takeLog10OfOutput ) x = Math.log10(x);
            normalized[i] = x;
        }

        return normalized;
    }
    /**
     * normalizes the log10-based array.  ASSUMES THAT ALL ARRAY ENTRIES ARE <= 0 (<= 1 IN REAL-SPACE).
     *
     * @param array  the array to be normalized
     *
     * @return a newly allocated array corresponding the normalized values in array
    */
    public static double[] normalizeFromLog10(double[] array) {
        return normalizeFromLog10(array, false);
    }

    public static double[] normalizeFromLog10(List<Double> array) {
        return normalizeFromLog10(array, false);
    }

    public static int maxElementIndex(double[] array) {
        if ( array == null ) throw new IllegalArgumentException("Array cannot be null!");

        int maxI = -1;
        for ( int i = 0; i < array.length; i++ ) {
            if ( maxI == -1 || array[i] > array[maxI] )
                maxI = i;
        }

        return maxI;
    }

    public static double arrayMax(double[] array) {
        return array[maxElementIndex(array)];
    }

    public static double arrayMin(double[] array) {
        return array[minElementIndex(array)];
    }

    public static byte arrayMin(byte[] array) {
        return array[minElementIndex(array)];
    }

    public static int minElementIndex(double[] array) {
        if ( array == null ) throw new IllegalArgumentException("Array cannot be null!");

        int minI = -1;
        for ( int i = 0; i < array.length; i++ ) {
            if ( minI == -1 || array[i] < array[minI] )
                minI = i;
        }

        return minI;
    }

    public static int minElementIndex(byte[] array) {
        if ( array == null ) throw new IllegalArgumentException("Array cannot be null!");

        int minI = -1;
        for ( int i = 0; i < array.length; i++ ) {
            if ( minI == -1 || array[i] < array[minI] )
                minI = i;
        }

        return minI;
    }    

    public static int arrayMax(List<Integer> array) {
        if ( array == null ) throw new IllegalArgumentException("Array cannot be null!");
        if ( array.size() == 0 ) throw new IllegalArgumentException("Array size cannot be 0!");

        int m = array.get(0);
        for ( int e : array ) m = Math.max(m, e);
        return m;
    }

    public static double arrayMax(List<Double> array) {
        if ( array == null ) throw new IllegalArgumentException("Array cannot be null!");
        if ( array.size() == 0 ) throw new IllegalArgumentException("Array size cannot be 0!");

        double m = array.get(0);
        for ( double e : array ) m = Math.max(m, e);
        return m;
    }

    public static double average(List<Long> vals, int maxI) {
        long sum = 0L;

        int i = 0;
        for (long x : vals) {
            if (i > maxI)
                break;
            sum += x;
            i++;
            //System.out.printf(" %d/%d", sum, i);
        }

        //System.out.printf("Sum = %d, n = %d, maxI = %d, avg = %f%n", sum, i, maxI, (1.0 * sum) / i);

        return (1.0 * sum) / i;
    }

    public static double averageDouble(List<Double> vals, int maxI) {
        double sum = 0.0;

        int i = 0;
        for (double x : vals) {
            if (i > maxI)
                break;
            sum += x;
            i++;
        }
        return (1.0 * sum) / i;
    }

    public static double average(List<Long> vals) {
        return average(vals, vals.size());
    }

    public static double averageDouble(List<Double> vals) {
        return averageDouble(vals, vals.size());
    }

    // Java Generics can't do primitive types, so I had to do this the simplistic way

    public static Integer[] sortPermutation(final int[] A) {
        class comparator implements Comparator<Integer> {
            public int compare(Integer a, Integer b) {
                if (A[a.intValue()] < A[b.intValue()]) {
                    return -1;
                }
                if (A[a.intValue()] == A[b.intValue()]) {
                    return 0;
                }
                if (A[a.intValue()] > A[b.intValue()]) {
                    return 1;
                }
                return 0;
            }
        }
        Integer[] permutation = new Integer[A.length];
        for (int i = 0; i < A.length; i++) {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }

    public static Integer[] sortPermutation(final double[] A) {
        class comparator implements Comparator<Integer> {
            public int compare(Integer a, Integer b) {
                if (A[a.intValue()] < A[b.intValue()]) {
                    return -1;
                }
                if (A[a.intValue()] == A[b.intValue()]) {
                    return 0;
                }
                if (A[a.intValue()] > A[b.intValue()]) {
                    return 1;
                }
                return 0;
            }
        }
        Integer[] permutation = new Integer[A.length];
        for (int i = 0; i < A.length; i++) {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }

    public static <T extends Comparable> Integer[] sortPermutation(List<T> A) {
        final Object[] data = A.toArray();

        class comparator implements Comparator<Integer> {
            public int compare(Integer a, Integer b) {
                return ((T) data[a]).compareTo(data[b]);
            }
        }
        Integer[] permutation = new Integer[A.size()];
        for (int i = 0; i < A.size(); i++) {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }


    public static int[] permuteArray(int[] array, Integer[] permutation) {
        int[] output = new int[array.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static double[] permuteArray(double[] array, Integer[] permutation) {
        double[] output = new double[array.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static Object[] permuteArray(Object[] array, Integer[] permutation) {
        Object[] output = new Object[array.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static String[] permuteArray(String[] array, Integer[] permutation) {
        String[] output = new String[array.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static <T> List<T> permuteList(List<T> list, Integer[] permutation) {
        List<T> output = new ArrayList<T>();
        for (int i = 0; i < permutation.length; i++) {
            output.add(list.get(permutation[i]));
        }
        return output;
    }


    /** Draw N random elements from list. */
    public static <T> List<T> randomSubset(List<T> list, int N) {
        if (list.size() <= N) {
            return list;
        }

        int idx[] = new int[list.size()];
        for (int i = 0; i < list.size(); i++) {
            idx[i] = GenomeAnalysisEngine.getRandomGenerator().nextInt();
        }

        Integer[] perm = sortPermutation(idx);

        List<T> ans = new ArrayList<T>();
        for (int i = 0; i < N; i++) {
            ans.add(list.get(perm[i]));
        }

        return ans;
    }

    // lifted from the internet
    // http://www.cs.princeton.edu/introcs/91float/Gamma.java.html
    public static double logGamma(double x) {
        double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
        double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1)
                + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
                + 0.00120858003 / (x + 4) - 0.00000536382 / (x + 5);
        return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
    }

    public static double percentage(double x, double base) {
        return (base > 0 ? (x / base) * 100.0 : 0);
    }

    public static double percentage(int x, int base) {
        return (base > 0 ? ((double) x / (double) base) * 100.0 : 0);
    }

    public static double percentage(long x, long base) {
        return (base > 0 ? ((double) x / (double) base) * 100.0 : 0);
    }

    public static int countOccurrences(char c, String s) {
        int count = 0;
        for (int i = 0; i < s.length(); i++) {
            count += s.charAt(i) == c ? 1 : 0;
        }
        return count;
    }

    public static <T> int countOccurrences(T x, List<T> l) {
        int count = 0;
        for (T y : l) {
            if (x.equals(y)) count++;
        }

        return count;
    }

    /**
     * Returns n random indices drawn with replacement from the range 0..(k-1)
     *
     * @param n the total number of indices sampled from
     * @param k the number of random indices to draw (with replacement)
     * @return a list of k random indices ranging from 0 to (n-1) with possible duplicates
     */
    static public ArrayList<Integer> sampleIndicesWithReplacement(int n, int k) {

        ArrayList<Integer> chosen_balls = new ArrayList <Integer>(k);
        for (int i=0; i< k; i++) {
            //Integer chosen_ball = balls[rand.nextInt(k)];
            chosen_balls.add(GenomeAnalysisEngine.getRandomGenerator().nextInt(n));
            //balls.remove(chosen_ball);
        }

        return chosen_balls;
    }

    /**
     * Returns n random indices drawn without replacement from the range 0..(k-1)
     *
     * @param n the total number of indices sampled from
     * @param k the number of random indices to draw (without replacement)
     * @return a list of k random indices ranging from 0 to (n-1) without duplicates
     */
    static public ArrayList<Integer> sampleIndicesWithoutReplacement(int n, int k) {
        ArrayList<Integer> chosen_balls = new ArrayList<Integer>(k);

        for (int i = 0; i < n; i++) {
            chosen_balls.add(i);
        }

        Collections.shuffle(chosen_balls, GenomeAnalysisEngine.getRandomGenerator());

        //return (ArrayList<Integer>) chosen_balls.subList(0, k);
        return new ArrayList<Integer>(chosen_balls.subList(0, k));
    }

    /**
     * Given a list of indices into a list, return those elements of the list with the possibility of drawing list elements multiple times

     * @param indices  the list of indices for elements to extract
     * @param list     the list from which the elements should be extracted
     * @param <T>      the template type of the ArrayList
     * @return         a new ArrayList consisting of the elements at the specified indices
     */
    static public <T> ArrayList<T> sliceListByIndices(List<Integer> indices, List<T> list) {
        ArrayList<T> subset = new ArrayList<T>();

        for (int i : indices) {
            subset.add(list.get(i));
        }

        return subset;
    }

    public static Comparable orderStatisticSearch(int orderStat, List<Comparable> list) {
        // this finds the order statistic of the list (kth largest element)
        // the list is assumed *not* to be sorted

        final Comparable x = list.get(orderStat);
        ListIterator iterator = list.listIterator();
        ArrayList lessThanX = new ArrayList();
        ArrayList equalToX = new ArrayList();
        ArrayList greaterThanX = new ArrayList();

        for(Comparable y : list) {
            if(x.compareTo(y) > 0) {
                lessThanX.add(y);
            } else if(x.compareTo(y) < 0) {
                greaterThanX.add(y);
            } else
                equalToX.add(y);
        }

        if(lessThanX.size() > orderStat)
            return orderStatisticSearch(orderStat, lessThanX);
        else if(lessThanX.size() + equalToX.size() >= orderStat)
            return orderStat;
        else
            return orderStatisticSearch(orderStat - lessThanX.size() - equalToX.size(), greaterThanX);

    }


    public static Object getMedian(List<Comparable> list) {
        return orderStatisticSearch((int) Math.ceil(list.size()/2), list);
    }

    public static byte getQScoreOrderStatistic(List<SAMRecord> reads, List<Integer> offsets, int k) {
        // version of the order statistic calculator for SAMRecord/Integer lists, where the
        // list index maps to a q-score only through the offset index
        // returns the kth-largest q-score.

        if( reads.size() == 0) {
            return 0;
        }

        ArrayList lessThanQReads = new ArrayList();
        ArrayList equalToQReads = new ArrayList();
        ArrayList greaterThanQReads = new ArrayList();
        ArrayList lessThanQOffsets = new ArrayList();
        ArrayList greaterThanQOffsets = new ArrayList();

        final byte qk = reads.get(k).getBaseQualities()[offsets.get(k)];

        for(int iter = 0; iter < reads.size(); iter ++) {
            SAMRecord read = reads.get(iter);
            int offset = offsets.get(iter);
            byte quality = read.getBaseQualities()[offset];

            if(quality < qk) {
                lessThanQReads.add(read);
                lessThanQOffsets.add(offset);
            } else if(quality > qk) {
                greaterThanQReads.add(read);
                greaterThanQOffsets.add(offset);
            } else {
                equalToQReads.add(reads.get(iter));
            }
        }

        if(lessThanQReads.size() > k)
            return getQScoreOrderStatistic(lessThanQReads, lessThanQOffsets, k);
        else if(equalToQReads.size() + lessThanQReads.size() >= k)
            return qk;
        else
            return getQScoreOrderStatistic(greaterThanQReads, greaterThanQOffsets, k - lessThanQReads.size() - equalToQReads.size());

    }

    public static byte getQScoreMedian(List<SAMRecord> reads, List<Integer> offsets) {
        return getQScoreOrderStatistic(reads, offsets, (int)Math.floor(reads.size()/2.));
    }

    // from http://en.wikipedia.org/wiki/Digamma_function
    // According to J.M. Bernardo AS 103 algorithm the digamma function for x, a real number, can be approximated by:
    public static double diGamma(final double x) {
        return Math.log(x) - ( 1.0 / (2.0 * x) )
                           - ( 1.0 / (12.0 * Math.pow(x, 2.0)) )
                           + ( 1.0 / (120.0 * Math.pow(x, 4.0)) )
                           - ( 1.0 / (252.0 * Math.pow(x, 6.0)) );
    }

    /** A utility class that computes on the fly average and standard deviation for a stream of numbers.
     * The number of observations does not have to be known in advance, and can be also very big (so that
     * it could overflow any naive summation-based scheme or cause loss of precision).
     * Instead, adding a new number <code>observed</code>
     * to a sample with <code>add(observed)</code> immediately updates the instance of this object so that
     * it contains correct mean and standard deviation for all the numbers seen so far. Source: Knuth, vol.2
     * (see also e.g. http://www.johndcook.com/standard_deviation.html for online reference).
     */
    public static class RunningAverage {
        private double mean = 0.0;
        private double s = 0.0;
        private long obs_count = 0;

        public void add(double obs) {
            obs_count++;
            double oldMean = mean;
            mean += ( obs - mean ) / obs_count; // update mean
            s += ( obs - oldMean ) * ( obs - mean );
        }

        public void addAll(Collection<Number> col) {
            for ( Number o : col ) {
                add(o.doubleValue());
            }
        }

        public double mean() { return mean; }
        public double stddev() { return Math.sqrt(s/(obs_count - 1)); }
        public double var() { return s/(obs_count - 1); }
        public long observationCount() { return obs_count; }

        public RunningAverage clone() {
            RunningAverage ra = new RunningAverage();
            ra.mean = this.mean;
            ra.s = this.s;
            ra.obs_count = this.obs_count;
            return ra;
        }

        public void merge(RunningAverage other) {
            if ( this.obs_count > 0 || other.obs_count > 0 ) { // if we have any observations at all
                this.mean = ( this.mean * this.obs_count + other.mean * other.obs_count ) / ( this.obs_count + other.obs_count );
                this.s += other.s;
            }
            this.obs_count += other.obs_count;
        }
    }
    
    //
    // useful common utility routines
    //
    public static double rate(long n, long d)           { return n / (1.0 * Math.max(d, 1)); }
    public static double rate(int n, int d)             { return n / (1.0 * Math.max(d, 1)); }

    public static long inverseRate(long n, long d)      { return n == 0 ? 0 : d / Math.max(n, 1); }
    public static long inverseRate(int n, int d)        { return n == 0 ? 0 : d / Math.max(n, 1); }

    public static double ratio(int num, int denom)      { return ((double)num) / (Math.max(denom, 1)); }
    public static double ratio(long num, long denom)    { return ((double)num) / (Math.max(denom, 1)); }

    public static final double[] log10Cache;
    public static final double[] jacobianLogTable;
    public static final int JACOBIAN_LOG_TABLE_SIZE = 101;
    public static final double JACOBIAN_LOG_TABLE_STEP = 0.1;
    public static final double INV_JACOBIAN_LOG_TABLE_STEP = 1.0/JACOBIAN_LOG_TABLE_STEP;
    public static final double MAX_JACOBIAN_TOLERANCE = 10.0;
    private static final int MAXN = 10000;

    static {
        log10Cache = new double[2*MAXN];
	    jacobianLogTable = new double[JACOBIAN_LOG_TABLE_SIZE];

        log10Cache[0] = Double.NEGATIVE_INFINITY;
        for (int k=1; k < 2*MAXN; k++)
            log10Cache[k] = Math.log10(k);

    	for (int k=0; k < JACOBIAN_LOG_TABLE_SIZE; k++) {
	        jacobianLogTable[k] = Math.log10(1.0+Math.pow(10.0,-((double)k)
						       * JACOBIAN_LOG_TABLE_STEP));

	    }
    }

    static public double softMax(final double[] vec) {
        double acc = vec[0];
        for (int k=1; k < vec.length; k++)
            acc = softMax(acc,vec[k]);

        return acc;

    }

    static public double max(double x0, double x1, double x2) {
        double a = Math.max(x0,x1);
        return Math.max(a,x2);
    }
    
    static public double softMax(final double x0, final double x1, final double x2) {
         // compute naively log10(10^x[0] + 10^x[1]+...)
         //        return Math.log10(MathUtils.sumLog10(vec));

         // better approximation: do Jacobian logarithm function on data pairs
         double a = softMax(x0,x1);
         return softMax(a,x2);
    }

     static public double softMax(final double x, final double y) {
         if (Double.isInfinite(x))
             return y;

         if (Double.isInfinite(y))
             return x;

         if (y >= x + MAX_JACOBIAN_TOLERANCE)
             return y;
         if (x >= y + MAX_JACOBIAN_TOLERANCE)
             return x;

         // OK, so |y-x| < tol: we use the following identity then:
         // we need to compute log10(10^x + 10^y)
         // By Jacobian logarithm identity, this is equal to
         // max(x,y) + log10(1+10^-abs(x-y))
         // we compute the second term as a table lookup
         // with integer quantization

         //double diff = Math.abs(x-y);
         double diff = x-y;
         double t1 =x;
         if (diff<0) { //
             t1 = y;
             diff= -diff;
         }
         // t has max(x,y), diff has abs(x-y)
         // we have pre-stored correction for 0,0.1,0.2,... 10.0
         //int ind = (int)Math.round(diff*INV_JACOBIAN_LOG_TABLE_STEP);
         int ind = (int)(diff*INV_JACOBIAN_LOG_TABLE_STEP+0.5);
         // gdebug+
         //double z =Math.log10(1+Math.pow(10.0,-diff));
         //System.out.format("x: %f, y:%f, app: %f, true: %f ind:%d\n",x,y,t2,z,ind);
         //gdebug-
         return t1+jacobianLogTable[ind];
         // return Math.log10(Math.pow(10.0,x) + Math.pow(10.0,y));
     }

}
