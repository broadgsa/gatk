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

import com.google.java.contract.Requires;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.math.BigDecimal;
import java.util.*;

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

    @Requires({"d > 0.0"})
    public static int fastPositiveRound(double d) {
        return (int) (d + 0.5);
    }

    public static int fastRound(double d) {
        if ( d > 0.0 ) {
            return fastPositiveRound(d);
        } else {
            return -1*fastPositiveRound(-1*d);
        }
    }

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
     * Calculates the log10 cumulative sum of an array with log10 probabilities
     * @param log10p the array with log10 probabilites
     * @param upTo index in the array to calculate the cumsum up to
     * @return the log10 of the cumulative sum
     */
    public static double log10CumulativeSumLog10(double [] log10p, int upTo) {
        return log10sumLog10(log10p, 0, upTo);
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
        return log10sumLog10(log10p, start, log10p.length);
    }

    public static double log10sumLog10(double[] log10p, int start, int finish) {
        double sum = 0.0;

        double maxValue = Utils.findMaxEntry(log10p);
        for ( int i = start; i < finish; i++ ) {
            sum += Math.pow(10.0, log10p[i] - maxValue);
        }

        return Math.log10(sum) + maxValue;
    }

    public static double sumDoubles(List<Double> values) {
        double s = 0.0;
        for ( double v : values) s += v;
        return s;
    }

    public static int sumIntegers(List<Integer> values) {
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

    public static double binomialCoefficient (int n, int k) {
        return Math.pow(10, log10BinomialCoefficient(n, k));
    }
    /**
     * Computes a binomial probability.  This is computed using the formula
     *
     *  B(k; n; p) = [ n! / ( k! (n - k)! ) ] (p^k)( (1-p)^k )
     *
     * where n is the number of trials, k is the number of successes, and p is the probability of success
     *
     * @param n  number of Bernoulli trials
     * @param k  number of successes
     * @param p  probability of success
     *
     * @return   the binomial probability of the specified configuration.  Computes values down to about 1e-237.
     */
    public static double binomialProbability (int n, int k, double p) {
        return Math.pow(10, log10BinomialProbability(n, k, Math.log10(p)));
    }

    /**
     * Performs the cumulative sum of binomial probabilities, where the probability calculation is done in log space.
     * @param start - start of the cumulant sum (over hits)
     * @param end - end of the cumulant sum (over hits)
     * @param total - number of attempts for the number of hits
     * @param probHit - probability of a successful hit
     * @return - returns the cumulative probability
     */
    public static double binomialCumulativeProbability(int start, int end, int total, double probHit) {
        double cumProb = 0.0;
        double prevProb;
        BigDecimal probCache = BigDecimal.ZERO;

        for(int hits = start; hits < end; hits++) {
            prevProb = cumProb;
            double probability = binomialProbability(total, hits, probHit);
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
     * Computes a multinomial coefficient efficiently avoiding overflow even for large numbers.
     * This is computed using the formula:
     *
     *   M(x1,x2,...,xk; n) = [ n! / (x1! x2! ... xk!) ]
     *
     * where xi represents the number of times outcome i was observed, n is the number of total observations.
     * In this implementation, the value of n is inferred as the sum over i of xi.
     *
     * @param k  an int[] of counts, where each element represents the number of times a certain outcome was observed
     * @return   the multinomial of the specified configuration.
     */
    public static double multinomialCoefficient (int [] k) {
        int n = 0;
        for (int xi : k) {
            n += xi;
        }

        return Math.pow(10, log10MultinomialCoefficient(n, k));
    }

    /**
     * Computes a multinomial probability efficiently avoiding overflow even for large numbers.
     * This is computed using the formula:
     *
     *   M(x1,x2,...,xk; n; p1,p2,...,pk) = [ n! / (x1! x2! ... xk!) ] (p1^x1)(p2^x2)(...)(pk^xk)
     *
     * where xi represents the number of times outcome i was observed, n is the number of total observations, and
     * pi represents the probability of the i-th outcome to occur.  In this implementation, the value of n is
     * inferred as the sum over i of xi.
     *
     * @param k  an int[] of counts, where each element represents the number of times a certain outcome was observed
     * @param p  a double[] of probabilities, where each element represents the probability a given outcome can occur
     * @return   the multinomial probability of the specified configuration.
     */
    public static double multinomialProbability (int[] k, double[] p) {
        if (p.length != k.length)
            throw new UserException.BadArgumentValue("p and k", "Array of log10 probabilities must have the same size as the array of number of sucesses: " + p.length + ", " + k.length);

        int n = 0;
        double [] log10P = new double[p.length];
        for (int i=0; i<p.length; i++) {
            log10P[i] = Math.log10(p[i]);
            n += k[i];
        }
        return Math.pow(10,log10MultinomialProbability(n, k, log10P));
    }

    /**
     * calculate the Root Mean Square of an array of integers
     * @param x  an byte[] of numbers
     * @return   the RMS of the specified numbers.
    */
    public static double rms(byte[] x) {
        if ( x.length == 0 )
            return 0.0;

        double rms = 0.0;
        for (int i : x)
            rms += i * i;
        rms /= x.length;
        return Math.sqrt(rms);
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

    public static double rms(Collection<Integer> l) {
        if (l.size() == 0)
            return 0.0;

        double rms = 0.0;
        for (int i : l)
            rms += i*i;
        rms /= l.size();
        return Math.sqrt(rms);
    }

    public static double distanceSquared( final double[] x, final double[] y ) {
        double dist = 0.0;
        for(int iii = 0; iii < x.length; iii++) {
            dist += (x[iii] - y[iii]) * (x[iii] - y[iii]);
        }
        return dist;
    }

    public static double round(double num, int digits) {
        double result = num * Math.pow(10.0, (double)digits);
        result = Math.round(result);
        result = result / Math.pow(10.0, (double)digits);
        return result;
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
        return normalizeFromLog10(array, takeLog10OfOutput, false);
    }

    public static double[] normalizeFromLog10(double[] array, boolean takeLog10OfOutput, boolean keepInLogSpace) {

        // for precision purposes, we need to add (or really subtract, since they're
        // all negative) the largest value; also, we need to convert to normal-space.
        double maxValue = Utils.findMaxEntry(array);

        // we may decide to just normalize in log space with converting to linear space
        if (keepInLogSpace) {
            for (int i = 0; i < array.length; i++)
                array[i] -= maxValue;
            return array;
        }

        // default case: go to linear space
        double[] normalized = new double[array.length];

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
        double maxValue = MathUtils.arrayMaxDouble( array );
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

    public static int maxElementIndex(int[] array) {
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

    public static int arrayMaxInt(List<Integer> array) {
        if ( array == null ) throw new IllegalArgumentException("Array cannot be null!");
        if ( array.size() == 0 ) throw new IllegalArgumentException("Array size cannot be 0!");

        int m = array.get(0);
        for ( int e : array ) m = Math.max(m, e);
        return m;
    }

    public static double arrayMaxDouble(List<Double> array) {
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

    public static byte average(byte[] vals) {
        int sum = 0;
        for (byte v : vals) {
            sum += v;
        }
        return (byte) Math.floor(sum/vals.length);
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

    public static int countOccurrences(byte element, byte [] array) {
        int count = 0;
        for (byte y : array) {
            if (element == y)
                count++;
        }

        return count;
    }

    /**
     * Returns the top (larger) N elements of the array. Naive n^2 implementation (Selection Sort).
     * Better than sorting if N (number of elements to return) is small
     *
     * @param array the array
     * @param n number of top elements to return
     * @return the n larger elements of the array
     */
    public static Collection<Double> getNMaxElements(double [] array, int n) {
        ArrayList<Double> maxN = new ArrayList<Double>(n);
        double lastMax = Double.MAX_VALUE;
        for (int i=0; i<n; i++) {
            double max = Double.MIN_VALUE;
            for (double x : array) {
                max = Math.min(lastMax, Math.max(x, max));
            }
            maxN.add(max);
            lastMax = max;
        }
        return maxN;
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
        // we need to compute log10(10^x + 10^y)
        // By Jacobian logarithm identity, this is equal to
        // max(x,y) + log10(1+10^-abs(x-y))
        // we compute the second term as a table lookup
        // with integer quantization

        // slow exact version:
        // return Math.log10(Math.pow(10.0,x) + Math.pow(10.0,y));

        double diff = x-y;

        if (diff > MAX_JACOBIAN_TOLERANCE)
            return x;
        else if (diff < -MAX_JACOBIAN_TOLERANCE)
            return y;
        else if (diff >= 0) {
            int ind = (int)(diff*INV_JACOBIAN_LOG_TABLE_STEP+0.5);
            return x + jacobianLogTable[ind];
        }
        else {
            int ind = (int)(-diff*INV_JACOBIAN_LOG_TABLE_STEP+0.5);
            return y + jacobianLogTable[ind];
        }
    }

    public static double phredScaleToProbability (byte q) {
        return Math.pow(10,(-q)/10.0);
    }

    public static double phredScaleToLog10Probability (byte q) {
        return ((-q)/10.0);
    }

    /**
     * Returns the phred scaled value of probability p
     * @param p probability (between 0 and 1).
     * @return phred scaled probability of p
     */
    public static byte probabilityToPhredScale (double p) {
        return (byte) ((-10) * Math.log10(p));
    }

    public static double log10ProbabilityToPhredScale (double log10p) {
        return (-10) * log10p;
    }

    /**
     * Converts LN to LOG10
     * @param ln log(x)
     * @return log10(x)
     */
    public static double lnToLog10 (double ln) {
        return ln * Math.log10(Math.exp(1));
    }

    /**
     * Constants to simplify the log gamma function calculation.
     */
    private static final double
        zero = 0.0,
        one  = 1.0,
        half = .5,
        a0  =  7.72156649015328655494e-02,
        a1  =  3.22467033424113591611e-01,
        a2  =  6.73523010531292681824e-02,
        a3  =  2.05808084325167332806e-02,
        a4  =  7.38555086081402883957e-03,
        a5  =  2.89051383673415629091e-03,
        a6  =  1.19270763183362067845e-03,
        a7  =  5.10069792153511336608e-04,
        a8  =  2.20862790713908385557e-04,
        a9  =  1.08011567247583939954e-04,
        a10 =  2.52144565451257326939e-05,
        a11 =  4.48640949618915160150e-05,
        tc  =  1.46163214496836224576e+00,
        tf  = -1.21486290535849611461e-01,
        tt  = -3.63867699703950536541e-18,
        t0  =  4.83836122723810047042e-01,
        t1  = -1.47587722994593911752e-01,
        t2  =  6.46249402391333854778e-02,
        t3  = -3.27885410759859649565e-02,
        t4  =  1.79706750811820387126e-02,
        t5  = -1.03142241298341437450e-02,
        t6  =  6.10053870246291332635e-03,
        t7  = -3.68452016781138256760e-03,
        t8  =  2.25964780900612472250e-03,
        t9  = -1.40346469989232843813e-03,
        t10 =  8.81081882437654011382e-04,
        t11 = -5.38595305356740546715e-04,
        t12 =  3.15632070903625950361e-04,
        t13 = -3.12754168375120860518e-04,
        t14 =  3.35529192635519073543e-04,
        u0  = -7.72156649015328655494e-02,
        u1  =  6.32827064025093366517e-01,
        u2  =  1.45492250137234768737e+00,
        u3  =  9.77717527963372745603e-01,
        u4  =  2.28963728064692451092e-01,
        u5  =  1.33810918536787660377e-02,
        v1  =  2.45597793713041134822e+00,
        v2  =  2.12848976379893395361e+00,
        v3  =  7.69285150456672783825e-01,
        v4  =  1.04222645593369134254e-01,
        v5  =  3.21709242282423911810e-03,
        s0  = -7.72156649015328655494e-02,
        s1  =  2.14982415960608852501e-01,
        s2  =  3.25778796408930981787e-01,
        s3  =  1.46350472652464452805e-01,
        s4  =  2.66422703033638609560e-02,
        s5  =  1.84028451407337715652e-03,
        s6  =  3.19475326584100867617e-05,
        r1  =  1.39200533467621045958e+00,
        r2  =  7.21935547567138069525e-01,
        r3  =  1.71933865632803078993e-01,
        r4  =  1.86459191715652901344e-02,
        r5  =  7.77942496381893596434e-04,
        r6  =  7.32668430744625636189e-06,
        w0  =  4.18938533204672725052e-01,
        w1  =  8.33333333333329678849e-02,
        w2  = -2.77777777728775536470e-03,
        w3  =  7.93650558643019558500e-04,
        w4  = -5.95187557450339963135e-04,
        w5  =  8.36339918996282139126e-04,
        w6  = -1.63092934096575273989e-03;

    /**
     * Efficient rounding functions to simplify the log gamma function calculation
     *   double to long with 32 bit shift
     */
    private static final int HI (double x) {
        return (int)(Double.doubleToLongBits(x) >> 32);
    }

    /**
     * Efficient rounding functions to simplify the log gamma function calculation
     *   double to long without shift
     */
    private static final int LO (double x) {
        return (int)Double.doubleToLongBits(x);
    }

    /**
     * Most efficent implementation of the lnGamma (FDLIBM)
     * Use via the log10Gamma wrapper method.
     */
    private static double lnGamma (double x) {
        double t,y,z,p,p1,p2,p3,q,r,w;
        int i;

        int hx = HI(x);
        int lx = LO(x);

        /* purge off +-inf, NaN, +-0, and negative arguments */
        int ix = hx&0x7fffffff;
        if (ix >= 0x7ff00000) return Double.POSITIVE_INFINITY;
        if ((ix|lx)==0 || hx < 0) return Double.NaN;
        if (ix<0x3b900000) {	/* |x|<2**-70, return -log(|x|) */
            return -Math.log(x);
        }

        /* purge off 1 and 2 */
        if((((ix-0x3ff00000)|lx)==0)||(((ix-0x40000000)|lx)==0)) r = 0;
        /* for x < 2.0 */
        else if(ix<0x40000000) {
            if(ix<=0x3feccccc) { 	/* lgamma(x) = lgamma(x+1)-log(x) */
                r = -Math.log(x);
                if(ix>=0x3FE76944) {y = one-x; i= 0;}
                else if(ix>=0x3FCDA661) {y= x-(tc-one); i=1;}
                else {y = x; i=2;}
            } else {
                r = zero;
                if(ix>=0x3FFBB4C3) {y=2.0-x;i=0;} /* [1.7316,2] */
                else if(ix>=0x3FF3B4C4) {y=x-tc;i=1;} /* [1.23,1.73] */
                else {y=x-one;i=2;}
            }

            switch(i) {
            case 0:
                z = y*y;
                p1 = a0+z*(a2+z*(a4+z*(a6+z*(a8+z*a10))));
                p2 = z*(a1+z*(a3+z*(a5+z*(a7+z*(a9+z*a11)))));
                p  = y*p1+p2;
                r  += (p-0.5*y); break;
            case 1:
                z = y*y;
                w = z*y;
                p1 = t0+w*(t3+w*(t6+w*(t9 +w*t12)));	/* parallel comp */
                p2 = t1+w*(t4+w*(t7+w*(t10+w*t13)));
                p3 = t2+w*(t5+w*(t8+w*(t11+w*t14)));
                p  = z*p1-(tt-w*(p2+y*p3));
                r += (tf + p); break;
            case 2:
                p1 = y*(u0+y*(u1+y*(u2+y*(u3+y*(u4+y*u5)))));
                p2 = one+y*(v1+y*(v2+y*(v3+y*(v4+y*v5))));
                r += (-0.5*y + p1/p2);
            }
        }
        else if(ix<0x40200000) { 			/* x < 8.0 */
            i = (int)x;
            t = zero;
            y = x-(double)i;
            p = y*(s0+y*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))));
            q = one+y*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))));
            r = half*y+p/q;
            z = one;	/* lgamma(1+s) = log(s) + lgamma(s) */
            switch(i) {
            case 7: z *= (y+6.0);	/* FALLTHRU */
            case 6: z *= (y+5.0);	/* FALLTHRU */
            case 5: z *= (y+4.0);	/* FALLTHRU */
            case 4: z *= (y+3.0);	/* FALLTHRU */
            case 3: z *= (y+2.0);	/* FALLTHRU */
                r += Math.log(z); break;
            }
            /* 8.0 <= x < 2**58 */
        } else if (ix < 0x43900000) {
            t = Math.log(x);
            z = one/x;
            y = z*z;
            w = w0+z*(w1+y*(w2+y*(w3+y*(w4+y*(w5+y*w6)))));
            r = (x-half)*(t-one)+w;
        } else
            /* 2**58 <= x <= inf */
            r =  x*(Math.log(x)-one);
        return r;
    }

    /**
     * Calculates the log10 of the gamma function for x using the efficient FDLIBM
     * implementation to avoid overflows and guarantees high accuracy even for large
     * numbers.
     *
     * @param x the x parameter
     * @return the log10 of the gamma function at x.
     */
    public static double log10Gamma (double x) {
        return lnToLog10(lnGamma(x));
    }

    /**
     * Calculates the log10 of the binomial coefficient. Designed to prevent
     * overflows even with very large numbers.
     *
     * @param n total number of trials
     * @param k number of successes
     * @return the log10 of the binomial coefficient
     */
    public static double log10BinomialCoefficient (int n, int k) {
        return log10Gamma(n+1) - log10Gamma(k+1) - log10Gamma(n-k+1);
    }

    public static double log10BinomialProbability (int n, int k, double log10p) {
        double log10OneMinusP = Math.log10(1-Math.pow(10,log10p));
        return log10BinomialCoefficient(n, k) + log10p*k + log10OneMinusP*(n-k);
    }


    /**
     * Calculates the log10 of the multinomial coefficient. Designed to prevent
     * overflows even with very large numbers.
     *
     * @param n total number of trials
     * @param k array of any size with the number of successes for each grouping (k1, k2, k3, ..., km)
     * @return
     */
    public static double log10MultinomialCoefficient (int n, int [] k) {
        double denominator = 0.0;
        for (int x : k) {
            denominator += log10Gamma(x+1);
        }
        return log10Gamma(n+1) - denominator;
    }

    /**
     * Computes the log10 of the multinomial distribution probability given a vector
     * of log10 probabilities. Designed to prevent overflows even with very large numbers.
     *
     * @param n number of trials
     * @param k array of number of successes for each possibility
     * @param log10p array of log10 probabilities
     * @return
     */
    public static double log10MultinomialProbability (int n, int [] k, double [] log10p) {
        if (log10p.length != k.length)
            throw new UserException.BadArgumentValue("p and k", "Array of log10 probabilities must have the same size as the array of number of sucesses: " + log10p.length + ", " + k.length);
        double log10Prod = 0.0;
        for (int i=0; i<log10p.length; i++) {
            log10Prod += log10p[i]*k[i];
        }
        return log10MultinomialCoefficient(n, k) + log10Prod;
    }

    public static double factorial (int x) {
        return Math.pow(10, log10Gamma(x+1));
    }

    public static double log10Factorial (int x) {
        return log10Gamma(x+1);
    }
}
