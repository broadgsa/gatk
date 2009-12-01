package org.broadinstitute.sting.utils;

import java.util.*;

public class WilcoxonRankSum {

    private LinkedList<Pair<Double, WILCOXON_SET>> observations = new LinkedList<Pair<Double, WILCOXON_SET>>();

    public enum WILCOXON_SET { SET1, SET2 }

    public WilcoxonRankSum() {}

    public void addObservation(Double observation, WILCOXON_SET set) {
        observations.add(new Pair<Double, WILCOXON_SET>(observation, set));
    }

    public double getTwoTailedPValue() {
        if ( observations.size() == 0 )
            return 0.0;

////////
// Remove me
////////
        observations.clear();

        for (int i=0; i < 50; i++)
            addObservation(1.0, WILCOXON_SET.SET2);
        for (int i=0; i < 50; i++)
            addObservation(3.0, WILCOXON_SET.SET1);

        observations.clear();

        observations.add(new Pair<Double, WILCOXON_SET>(1.0, WILCOXON_SET.SET2));
        observations.add(new Pair<Double, WILCOXON_SET>(1.0, WILCOXON_SET.SET2));
        observations.add(new Pair<Double, WILCOXON_SET>(1.0, WILCOXON_SET.SET2));
        observations.add(new Pair<Double, WILCOXON_SET>(1.0, WILCOXON_SET.SET2));

        observations.add(new Pair<Double, WILCOXON_SET>(8.0, WILCOXON_SET.SET1));
        observations.add(new Pair<Double, WILCOXON_SET>(8.0, WILCOXON_SET.SET1));
        observations.add(new Pair<Double, WILCOXON_SET>(8.0, WILCOXON_SET.SET1));
        observations.add(new Pair<Double, WILCOXON_SET>(8.0, WILCOXON_SET.SET1));

        observations.clear();

        observations.add(new Pair<Double, WILCOXON_SET>(2.0, WILCOXON_SET.SET1));
        observations.add(new Pair<Double, WILCOXON_SET>(4.0, WILCOXON_SET.SET1));
        observations.add(new Pair<Double, WILCOXON_SET>(6.0, WILCOXON_SET.SET1));
        observations.add(new Pair<Double, WILCOXON_SET>(8.0, WILCOXON_SET.SET1));

        observations.add(new Pair<Double, WILCOXON_SET>(1.0, WILCOXON_SET.SET2));
        observations.add(new Pair<Double, WILCOXON_SET>(2.0, WILCOXON_SET.SET2));
        observations.add(new Pair<Double, WILCOXON_SET>(3.0, WILCOXON_SET.SET2));
        observations.add(new Pair<Double, WILCOXON_SET>(4.0, WILCOXON_SET.SET2));
////////


        // sort
        Collections.sort(observations, new PairComparator());

        // rank
        double[] ranks = calculateRanks(observations);
        // for (int i = 0; i < ranks.length; i++)
        //     System.out.println(observations.get(i).first + " -> " + ranks[i]);

        // sum
        double sum = 0.0;
        int n1 = 0;
        for (int i = 0; i < ranks.length; i++) {
            if ( observations.get(i).second == WILCOXON_SET.SET1 ) {
                sum += ranks[i];
                n1++;
            }
        }
        int n2 = ranks.length - n1;

        // we want the smaller of U1 and U2
        double U1 = sum - (n1 * (n1 + 1.0) / 2.0);
        double U2 = (n1 * n2) - U1;                                                                
        double U = Math.min(U1, U2);

        // calculate the normal approximation
        double MuU = n1 * n2 / 2.0;
        double stdevU = Math.sqrt(n1 * n2 * (n1 + n2 + 1.0) / 12.0);
        double z = (U - MuU) / stdevU;

        // System.out.println("U1=" + U1);
        // System.out.println("U2=" + U2);
        // System.out.println("U=" + U);
        // System.out.println("z=" + z);

        // compute p-value
        // double pvalue = ndtr(z);  // normal distribution function
        double pvalue = 0.0;

        return pvalue;
    }

    private static double[] calculateRanks(List<Pair<Double, WILCOXON_SET>> observations) {
        int length = observations.size();
        double[] ranks = new double[length];

        int currentIndex = 1;
        Double currentValue = observations.get(0).first;
        int startIndex = 0;
        int endIndex = 0;

        while ( currentIndex < length ) {
            // if two observations have the same value, they'll need to be ranked together
            if ( observations.get(currentIndex).first.equals(currentValue) ) {
                endIndex++;
            } else {
                setRanks(ranks, startIndex, endIndex);

                // increment the holders
                startIndex = endIndex = currentIndex;
                currentValue = observations.get(currentIndex).first;
            }
            currentIndex++;
        }
        if ( startIndex < length )
            setRanks(ranks, startIndex, endIndex);

        return ranks;
    }

    private static void setRanks(double[] ranks, int startIndex, int endIndex) {

        // the rank is the average rank of all equal observations
        double rankValue = 0.0;
        for (int i = startIndex; i <= endIndex; i++)
            rankValue += (i+1);
        rankValue /= (endIndex - startIndex + 1);

        // set the value
        for (int i = startIndex; i <= endIndex; i++)
            ranks[i] = rankValue;
    }

    private class PairComparator implements Comparator<Pair<Double, WILCOXON_SET>>{
        public int compare(Pair<Double, WILCOXON_SET> p1, Pair<Double, WILCOXON_SET> p2) {
            return (p1.first < p2.first ? -1 : (p1.first == p2.first ? 0 : 1));
        }
    }
}