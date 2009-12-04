package org.broadinstitute.sting.utils;

import cern.jet.random.Normal;

import java.util.*;

public class WilcoxonRankSum {

    // *****************************************************************************************//
    // The following 4 variables were copied from Tim Fennell's RankSumTest.java code in Picard //
    // *****************************************************************************************//

    // Constructs a normal distribution; actual values of mean and SD don't matter since it's
    // just used to convert a z-score into a cumulative probability
    private static final double NORMAL_MEAN = 100;
    private static final double NORMAL_SD   = 15;
    private static final Normal NORMAL = new Normal(NORMAL_MEAN, NORMAL_SD, null);

    // The minimum length for both data series (individually) in order to use a normal distribution
    // to calculate the Z-score and the p-value. If either series is shorter than this value then
    // we don't attempt to assign a p-value
    private static final int minimumNormalN = 10;

    // *****************************************************************************************//

    public enum WILCOXON_SET { SET1, SET2 }

    // random number generator for dithering
    private static final long RANDOM_SEED = 1252863494;
    private Random generator = new Random(RANDOM_SEED);

    // storage for observations
    private LinkedList<Pair<Double, WILCOXON_SET>> observations = new LinkedList<Pair<Double, WILCOXON_SET>>();


    public WilcoxonRankSum() {}

    // add an observation for a given set
    public void addObservation(Double observation, WILCOXON_SET set) {
        observations.add(new Pair<Double, WILCOXON_SET>(observation, set));
    }

    // calculate normal approximation of the p-value
    // returns -1 when unable to calculate it (too few data points)
    public double getTwoTailedPValue() {
        if ( observations.size() == 0 )
            return -1.0;

        // dither to break rank ties
        dither();

        // sort
        Collections.sort(observations, new PairComparator());

        // sum
        double sum = 0.0;
        int n1 = 0;
        for (int i = 0; i < observations.size(); i++) {
            if ( observations.get(i).second == WILCOXON_SET.SET1 ) {
                sum += i+1;
                n1++;
            }
        }
        int n2 = observations.size() - n1;

        // if we don't have enough data points, quit
        if ( n1 < minimumNormalN || n2 < minimumNormalN )
            return -1.0;
                
        // we want the smaller of U1 and U2
        double U1 = sum - (n1 * (n1 + 1.0) / 2.0);
        double U2 = (n1 * n2) - U1;                                                                
        double U = Math.min(U1, U2);

        // calculate the normal approximation
        double MuU = n1 * n2 / 2.0;
        double stdevU = Math.sqrt(n1 * n2 * (n1 + n2 + 1.0) / 12.0);
        double z = (U - MuU) / stdevU;

        // compute p-value.  Taken from Tim Fennell's RankSumTest.java code in Picard
        double pvalue = 2.0 * NORMAL.cdf(NORMAL_MEAN + z * NORMAL_SD);

        // for (int i = 0; i < observations.size(); i++)
        //      System.out.println(observations.get(i).first + " -> set" + (observations.get(i).second == WILCOXON_SET.SET1 ? 1 : 2));
        // System.out.println("U1=" + U1);
        // System.out.println("U2=" + U2);
        // System.out.println("U=" + U);
        // System.out.println("Zscore=" + z);
        // System.out.println("Pvalue=" + pvalue);

        return pvalue;
    }

    private void dither() {
        for ( Pair<Double, WILCOXON_SET> observation : observations ) {
            // generate a random number between 0 and 10,000
            int rand = generator.nextInt(10000);

            // convert it into a small floating point number by dividing by 1,000,000
            double smallFloat = (double)rand / 1000000;

            // add it to the observation
            observation.first += smallFloat;
        }
    }
    
    private class PairComparator implements Comparator<Pair<Double, WILCOXON_SET>>{
        public int compare(Pair<Double, WILCOXON_SET> p1, Pair<Double, WILCOXON_SET> p2) {
            return (p1.first < p2.first ? -1 : (p1.first == p2.first ? 0 : 1));
        }
    }
}