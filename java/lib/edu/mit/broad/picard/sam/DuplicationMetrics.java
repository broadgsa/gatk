package edu.mit.broad.picard.sam;

import edu.mit.broad.picard.metrics.MetricBase;
import edu.mit.broad.picard.util.Histogram;

/**
 * Metrics that are calculated during the process of marking duplicates
 * within a stream of SAMRecords.
 */
public class DuplicationMetrics extends MetricBase {
    /** The number of mapped reads examined which did not have a mapped mate pair. */
    public long UNPAIRED_READS_EXAMINED;

    /** The number of mapped read pairs examined. */
    public long READ_PAIRS_EXAMINED;

    /** The total number of unmapped reads examined. */
    public long UNMAPPED_READS;

    /** The number of fragments that were marked as duplicates. */
    public long UNPAIRED_READ_DUPLICATES;

    /** The number of read pairs that were marked as duplicates. */
    public long READ_PAIR_DUPLICATES;

    /** The percentage of mapped sequence that is marked as duplicate. */
    public Double PERCENT_DUPLICATION;

    /** The estimated number of unique molecules in the library based on PE duplication. */
    public Long ESTIMATED_LIBRARY_SIZE;

    /**
     * Fills in the ESTIMATED_LIBRARY_SIZE based on the paired read data examined where
     * possible and the PERCENT_DUPLICATION.
     */
    public void calculateDerivedMetrics() {
        if (READ_PAIRS_EXAMINED > 0) {
            // Following code "borrowed" from CRD codebase
            long n    = READ_PAIRS_EXAMINED;
            long c = READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES;

            double m = 1.0, M = 100.0;

            if (c >= n || f(m*c, c, n) <= 0) {
                throw new IllegalStateException("Invalid values for pairs and unique pairs: "
                        + n + ", " + c);

            }

            while( f(M*c, c, n) >= 0 ) M *= 10.0;

            for (int i=0; i<40; i++ ) {
                double r = (m+M)/2.0;
                double u = f( r * c, c, n );
                if ( u == 0 ) break;
                else if ( u > 0 ) m = r;
                else if ( u < 0 ) M = r;
            }

            this.ESTIMATED_LIBRARY_SIZE = (long) (c * (m+M)/2.0);
        }

        PERCENT_DUPLICATION = (UNPAIRED_READ_DUPLICATES + READ_PAIR_DUPLICATES *2) /(double) (UNPAIRED_READS_EXAMINED + READ_PAIRS_EXAMINED *2);
    }

    /** Method that is used in the computation of estimated library size. */
    private double f(double x, double c, double n) {
        return c/x - 1 + Math.exp(-n/x);
    }

    /**
     * Estimates the ROI (return on investment) that one would see if a library was sequenced to
     * x higher coverage than the observed coverage.
     *
     * @param estimatedLibrarySize the estimated number of molecules in the library
     * @param x the multiple of sequencing to be simulated (i.e. how many X sequencing)
     * @param pairs the number of pairs observed in the actual sequencing
     * @param uniquePairs the number of unique pairs observed in the actual sequencing
     * @return a number z <= x that estimates if you had pairs*x as your sequencing then you
     *         would observe uniquePairs*z unique pairs.
     */
    private double estimateRoi(long estimatedLibrarySize, double x, long pairs, long uniquePairs) {
        return estimatedLibrarySize * ( 1 - Math.exp(-(x*pairs)/estimatedLibrarySize) ) / uniquePairs;
    }

    /**
     * Calculates a histogram using the estimateRoi method to estimate the effective yield
     * doing x sequencing for x=1..10.
     */
    public Histogram<Double> calculateRoiHistogram() {
        if (ESTIMATED_LIBRARY_SIZE == null) {
            try { calculateDerivedMetrics();  }
            catch (IllegalStateException ise) { return null; }
        }

        long uniquePairs = READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES;
        Histogram<Double> histo = new Histogram<Double>();

        for (double x=1; x<=10; x+=1) {
            histo.increment(x, estimateRoi(ESTIMATED_LIBRARY_SIZE, x, READ_PAIRS_EXAMINED, uniquePairs));
        }

        return histo;
    }

    // Main method used for debugging the derived metrics
//    public static void main(String[] args) {
//        DuplicationMetrics m = new DuplicationMetrics();
//        m.PAIRS_EXAMINED  = Integer.parseInt(args[0]);
//        m.DUPLICATE_PAIRS = m.PAIRS_EXAMINED - Integer.parseInt(args[1]);
//        m.calculateDerivedMetrics();
//        System.out.println("Percent Duplication: " + m.PERCENT_DUPLICATION);
//        System.out.println("Est. Library Size  : " + m.ESTIMATED_LIBRARY_SIZE);
//        System.out.println(m.calculateRoiHistogram());
//    }
}
