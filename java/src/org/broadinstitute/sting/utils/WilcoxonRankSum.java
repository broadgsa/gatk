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

import cern.jet.random.Normal;

import java.util.*;

import org.broadinstitute.sting.utils.collections.Pair;

public class WilcoxonRankSum {
    static final String headerString = ("nA nB .005 .01 .025 .05 .10 .20");
    static final String header[] = headerString.split(" ");
    static final double sigs[] = {-1, -1, .005, .01, .025, .05, .10, .20};

    // Probabilities relate to the distribution of WA, the rank sum for group A when
    // H0 : A = B is true. The tabulated value for the lower tail is the largest
    // value of wA for which pr(WA <= wA) <= prob . The tabulated value for the
    // upper tail is the smallest value of wfor which pr(W >= w) <= prob .

    // I think this data is wrong -- we should really do the computation outselves and remember the results
    static final int data[][] = {
            {4,     4,      -1,     -1,     10,     11,     13,     14},
            {4,     5,      -1,     10,     11,     12,     14,     15},
            {4,     6,      10,     11,     12,     13,     15,     17},
            {4,     7,      10,     11,     13,     14,     16,     18},
            {4,     8,      11,     12,     14,     15,     17,     20},
            {4,     9,      11,     13,     14,     16,     19,     21},
            {4,     10,     12,     13,     15,     17,     20,     23},
            {4,     11,     12,     14,     16,     18,     21,     24},
            {4,     12,     13,     15,     17,     19,     22,     26},
            {5,     5,      5,      15,     16,     17,     19,     20},
            {5,     6,      16,     17,     18,     20,     22,     24},
            {5,     7,      16,     18,     20,     21,     23,     26},
            {5,     8,      17,     19,     21,     23,     25,     28},
            {5,     9,      18,     20,     22,     24,     27,     30},
            {5,     10,     19,     21,     23,     26,     28,     32},
            {5,     11,     20,     22,     24,     27,     30,     34},
            {5,     12,     21,     23,     26,     28,     32,     36},
            {6,     6,      23,     24,     26,     28,     30,     33},
            {6,     7,      24,     25,     27,     29,     32,     35},
            {6,     8,      25,     27,     29,     31,     34,     37},
            {6,     9,      26,     28,     31,     33,     36,     40},
            {6,     10,     27,     29,     32,     35,     38,     42},
            {6,     11,     28,     30,     34,     37,     40,     44},
            {6,     12,     30,     32,     35,     38,     42,     47},
            {7,     7,      32,     34,     36,     39,     41,     45},
            {7,     8,      34,     35,     38,     41,     44,     48},
            {7,     9,      35,     37,     40,     43,     46,     50},
            {7,     10,     37,     39,     42,     45,     49,     53},
            {7,     11,     38,     40,     44,     47,     51,     56},
            {7,     12,     40,     42,     46,     49,     54,     59},
            {8,     8,      43,     45,     49,     51,     55,     59},
            {8,     9,      45,     47,     51,     54,     58,     62},
            {8,     10,     47,     49,     53,     56,     60,     65},
            {8,     11,     49,     51,     55,     59,     63,     69},
            {8,     12,     51,     53,     58,     62,     66,     72},
            {9,     9,      56,     59,     62,     66,     70,     75},
            {9,     10,     58,     61,     65,     69,     73,     78},
            {9,     11,     61,     63,     68,     72,     76,     82},
            {9,     12,     63,     66,     71,     75,     80,     86}};

    // *****************************************************************************************//
    // The following 4 variables were copied from Tim Fennell's RankSumTest.java code in Picard //
    // *****************************************************************************************//

    // Constructs a normal distribution; actual values of mean and SD don't matter since it's
    // just used to convert a z-score into a cumulative probability
    public boolean DEBUG = false;

    private static final double NORMAL_MEAN = 0;
    private static final double NORMAL_SD   = 1;
    private static final Normal NORMAL = new Normal(NORMAL_MEAN, NORMAL_SD, null);

    // The minimum length for both data series (individually) in order to use a normal distribution
    // to calculate the Z-score and the p-value. If either series is shorter than this value then
    // we don't attempt to assign a p-value
    private static final int minimumNormalN = 5;

    // *****************************************************************************************//

    public enum WILCOXON_SET { SET1, SET2 }
    public enum WILCOXON_H0 { SET1_NE_SET2, SET1_LT_SET2, SET1_GT_SET2, SMALLER_SET_LT }

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
    public double getPValue(WILCOXON_H0 h0) {
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

        // todo -- these are actually integers
        // we want the smaller of U1 and U2
        double U1 = sum - (n1 * (n1 + 1.0) / 2.0);
        double U2 = (n1 * n2) - U1;

        double pvalue;
        // if we don't have enough data points, quit
        if ( n1 < minimumNormalN || n2 < minimumNormalN ) {
            pvalue = exactCalculation(h0, n1, n2, U1, U2);
        } else {
            pvalue = normalApproximation(h0, n1, n2, U1, U2);
        }

        if ( DEBUG && (n1 < minimumNormalN || n2 < minimumNormalN) ) {
            //for (int i = 0; i < observations.size(); i++)
            //    System.out.println(observations.get(i).first + " -> set" + (observations.get(i).second == WILCOXON_SET.SET1 ? "Alt" : "Ref"));
            //System.out.printf("n1 %d n2 %d U1 %f U2 %f pValue %f QPValue %f%n", n1, n2, U1, U2, pvalue, QualityUtils.phredScaleErrorRate(pvalue));
        }

        return pvalue;
    }

    private double exactCalculation(WILCOXON_H0 h0, int n1, int n2, double U1, double U2) {
        double U = 0;
        switch (h0) {
//            case SET1_NE_SET2: // two-tailed
//                double U = Math.min(U1, U2);
//                z = (U - MuU) / stdevU;
//                pvalue = 2.0 * NORMAL.cdf(z);
//                break;
//            case SET1_LT_SET2: // one-tailed, SET1 < SET2
//                z = (U1 - MuU) / stdevU;
//                pvalue = NORMAL.cdf(z);
//                break;
//            case SET1_GT_SET2: // one-tailed, SET1 < SET2
//                z = (U2 - MuU) / stdevU;
//                pvalue = NORMAL.cdf(z);
//                break;
            case SMALLER_SET_LT: // one-tailed that the smaller set is < the bigger set
                U = n1 < n2 ? U1 : U2;
                break;
            default:
                throw new StingException("Unexpected WILCOXON H0: " + h0);
        }

        // data is nA nB then
        double pvalue = 1;  // if we don't find anything, the pvalue is 1 [not significant]
        for ( int nAi = 0; nAi < data.length; nAi++ ) {
            if ( data[nAi][0] == n1 && data[nAi][1] == n2 ) {
                for ( int j = 2; j < data[nAi].length; j++ ) {
                    if ( U < data[nAi][j] ) {
                        // we are significant
                        pvalue = sigs[j];
                        //System.out.printf("Found %d %d %f a match at %d %d with %d %f%n", n1, n2, U, nAi, j, data[nAi][j], pvalue );
                        break;
                    }
                }
            }
        }

        return pvalue;
    }

    private double normalApproximation(WILCOXON_H0 h0, int n1, int n2, double U1, double U2) {
        // calculate the normal approximation
        double MuU = n1 * n2 / 2.0;
        double stdevU = Math.sqrt(n1 * n2 * (n1 + n2 + 1.0) / 12.0);

        double z, pvalue;
        switch (h0) {
            case SET1_NE_SET2: // two-tailed
                double U = Math.min(U1, U2);
                z = (U - MuU) / stdevU;
                pvalue = 2.0 * NORMAL.cdf(z);
                break;
            case SET1_LT_SET2: // one-tailed, SET1 < SET2
                z = (U1 - MuU) / stdevU;
                pvalue = NORMAL.cdf(z);
                break;
            case SET1_GT_SET2: // one-tailed, SET1 < SET2
                z = (U2 - MuU) / stdevU;
                pvalue = NORMAL.cdf(z);
                break;
            case SMALLER_SET_LT: // one-tailed that the smaller set is < the bigger set
                double smallerU = n1 < n2 ? U1 : U2;
                z = (smallerU - MuU) / stdevU;
                pvalue = NORMAL.cdf(z);
                break;
            default:
                throw new StingException("Unexpected WILCOXON H0: " + h0);
        }

        return pvalue;
    }

    private void printValues() {
        for ( Pair<Double, WILCOXON_SET> obs : observations ) {
            System.out.printf("%s%n", obs);
        }
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