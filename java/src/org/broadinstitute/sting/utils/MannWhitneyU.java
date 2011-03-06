package org.broadinstitute.sting.utils;

import cern.jet.math.Arithmetic;
import cern.jet.random.Normal;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.Comparator;
import java.util.TreeSet;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 */
public class MannWhitneyU {

    private static Normal STANDARD_NORMAL = new Normal(0.0,1.0,null);

    private TreeSet<Pair<Number,USet>> observations;
    private int sizeSet1;
    private int sizeSet2;

    public MannWhitneyU() {
        observations = new TreeSet<Pair<Number,USet>>(new DitheringComparator());
        sizeSet1 = 0;
        sizeSet2 = 0;
    }

    /**
     * Add an observation into the observation tree
     * @param n: the observation (a number)
     * @param set: whether the observation comes from set 1 or set 2
     */
    public void add(Number n, USet set) {
        observations.add(new Pair<Number,USet>(n,set));
        if ( set == USet.SET1 ) {
            ++sizeSet1;
        } else {
            ++sizeSet2;
        }
    }

    /**
     * temporary method that will be generalized. Runs the standard two-sided test,
     * returns the u and p values.
     * @Returns a pair holding the u and p-value.
     */
    public Pair<Integer,Double> runTwoSidedTest() {
        Pair<Integer,USet> uPair = calculateTwoSidedU(observations);
        int u = uPair.first;
        int n = uPair.second == USet.SET1 ? sizeSet1 : sizeSet2;
        int m = uPair.second == USet.SET1 ? sizeSet2 : sizeSet1;
        double pval = calculateP(n,m,u,true);
        return new Pair<Integer,Double>(u,pval);
    }

    /**
     * Given a u statistic, calculate the p-value associated with it, dispatching to approximations where appropriate
     * @param n - The number of entries in the DOMINATED set
     * @param m - The number of entries in the DOMINANT set
     * @param u - the Mann-Whitney U value
     * @param twoSided - is the test twosided
     * @return the (possibly approximate) p-value associated with the MWU test
     */
    public static double calculateP(int n, int m, int u, boolean twoSided) {
        double pval;
        if ( n > 8 && m > 8 ) {
            pval = calculatePNormalApproximation(n,m,u);
        } else if ( n > 4 && m > 7 ) {
            pval = calculatePUniformApproximation(n,m,u);
        } else {
            pval = calculatePRecursively(n,m,u);
        }

        return twoSided ? 2*pval : pval;
    }

    /**
     * Uses a normal approximation to the U statistic in order to return a cdf p-value. See Mann, Whitney [1947]
     * @param n - The number of entries in the DOMINATED set
     * @param m - The number of entries in the DOMINANT set
     * @param u - the Mann-Whitney U value
     * @return p-value associated with the normal approximation
     */
    public static double calculatePNormalApproximation(int n,int m,int u) {
        double mean = ((double) m*n+1)/2;
        double var = (n*m*(n+m+1))/12;
        double z = ( u - mean )/Math.sqrt(var);
        return z < 0 ? STANDARD_NORMAL.cdf(z) : 1.0-STANDARD_NORMAL.cdf(z);
    }

    /**
     * Uses a sum-of-uniform-0-1 random variable approximation to the U statistic in order to return an approximate
     * p-value. See Buckle, Kraft, van Eeden [1969] (approx) and Billingsly [1995] or Stephens [1966] (sum of uniform CDF)
     * @param n -
     * @param m -
     * @param u -
     * @return
     */
    public static double calculatePUniformApproximation(int n, int m, int u) {
        int R = u + (n*(n+1))/2;
        double a = Math.sqrt(m*(n+m+1));
        double b = (n/2.0)*(1-Math.sqrt((n+m+1)/m));
        double z = b + R/a;
        if ( z < 0 ) { return 0.0; }
        else if ( z > n ) { return 1.0; }
        else {
            return 1/((double)Arithmetic.factorial(n))*uniformSumHelper(z, (int) Math.floor(z), n, 0);
        }
    }

    /**
     * Helper function for the sum of n uniform random variables
     * @param z - value at which to compute the (un-normalized) cdf
     * @param m - a cutoff integer (defined by m <= z < m + 1)
     * @param n - the number of uniform random variables
     * @param k - holder variable for the recursion (alternatively, the index of the term in the sequence)
     * @return the (un-normalized) cdf for the sum of n random variables
     */
    private static double uniformSumHelper(double z, int m, int n, int k) {
        if ( k > m ) { return 0; }
        int coef = (k % 2 == 0) ? 1 : -1;
        return coef*Arithmetic.binomial(n,k)*Math.pow(z-k,n) + uniformSumHelper(z,m,n,k+1);
    }

    /**
     * Calculates the U-statistic associated with a two-sided test (e.g. the RV from which one set is drawn
     * stochastically dominates the RV from which the other set is drawn); two-sidedness is accounted for
     * later on simply by multiplying the p-value by 2
     * @param observed
     * @return the minimum of the U counts (set1 dominates 2, set 2 dominates 1)
     */
    public static Pair<Integer,USet> calculateTwoSidedU(TreeSet<Pair<Number,USet>> observed ) {
        int set1SeenSoFar = 0;
        int set2SeenSoFar = 0;
        int uSet1DomSet2 = 0;
        int uSet2DomSet1 = 0;
        USet previous = null;
        for ( Pair<Number,USet> dataPoint : observed ) {

            if ( dataPoint.second == USet.SET1 ) {
                ++set1SeenSoFar;
            } else {
                ++set2SeenSoFar;
            }

            if ( previous != null ) {
                if ( dataPoint.second == USet.SET1 ) {
                    uSet2DomSet1 += set2SeenSoFar;
                } else {
                    uSet1DomSet2 += set1SeenSoFar;
                }
            }

            previous = dataPoint.second;
        }

        return uSet1DomSet2 < uSet2DomSet1 ? new Pair<Integer,USet>(uSet1DomSet2,USet.SET1) : new Pair<Integer,USet>(uSet2DomSet1,USet.SET2);
    }

    /**
     * The Mann-Whitney U statistic follows a recursive equation (that enumerates the proportion of possible
     * binary strings of "n" zeros, and "m" ones, where a one precedes a zero "u" times). This accessor
     * calls into that recursive calculation.
     * @param n: number of set-one entries (hypothesis: set-one is dominated by set-two)
     * @param m: number of set-two entries
     * @param u: number of set-two entries that precede set-one entries (e.g. 0,1,0,1,0 -> 3 )
     * @return the probability under the hypothesis that all sequences are equally likely of finding a set-two entry preceding a set-one entry "u" times.
     */
    public static double calculatePRecursively(int n, int m, int u) {
        if ( m + n > 16 ) { throw new StingException("Please use the appropriate (normal or sum of uniform) approximation"); }
        return cpr(n,m,u);
    }

    /**
     * @doc: just a shorter name for calculatePRecursively. See Mann, Whitney, [1947]
     * @n: number of set-1 entries
     * @m: number of set-2 entries
     * @u: number of times a set-2 entry as preceded a set-1 entry
     */
    private static double cpr(int n, int m, int u) {
        if ( u < 0 || n == 0 && m == 0 ) {
            return 0.0;
        }
        if ( m == 0 || n == 0 ) {
            // there are entries in set 1 or set 2, so no set-2 entry can precede a set-1 entry; thus u must be zero.
            // note that this exists only for edification, as when we reach this point, the coefficient on this term is zero anyway
            return ( u == 0 ) ? 1.0 : 0.0;
        }


        return (((double)n)/(n+m))*cpr(n-1,m,u-m) + (((double)m)/(n+m))*cpr(n,m-1,u);
    }

    /**
     * A comparator class which uses dithering on tie-breaking to ensure that the internal treeset drops no values
     * and to ensure that rank ties are broken at random.
     */
    private class DitheringComparator implements Comparator<Pair<Number,USet>> {

        public DitheringComparator() {}

        public boolean equals(Object other) { return false; }

        public int compare(Pair<Number,USet> left, Pair<Number,USet> right) {
            double comp = Double.compare(left.first.doubleValue(),right.first.doubleValue());
            if ( comp > 0 ) { return 1; }
            if ( comp < 0 ) { return -1; }
            return MathUtils.rand.nextBoolean() ? -1 : 1;
        }
    }

    public enum USet { SET1, SET2 }

}
