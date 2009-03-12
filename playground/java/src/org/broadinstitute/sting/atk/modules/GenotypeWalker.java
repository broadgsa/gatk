package org.broadinstitute.sting.atk.modules;

import org.broadinstitute.sting.atk.LocusIterator;
import org.broadinstitute.sting.atk.GenotypeEvidence;
import org.broadinstitute.sting.atk.LocusContext;
import org.broadinstitute.sting.utils.ReferenceOrderedDatum;
import net.sf.samtools.SAMRecord;


import java.util.List;
import static java.lang.System.currentTimeMillis;

public class GenotypeWalker extends BasicLociWalker<Integer, Integer> {
    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        //char[] = new char(26);
        long start_tm = currentTimeMillis();
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        String bases = "";
        String quals = "";
            //String offsetString = "";
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

                //if ( offset >= read.getReadString().length() )
                //    System.out.printf("  [%2d] [%s] %s%n", offset, read.format(), read.getReadString());

            bases += read.getReadString().charAt(offset);
            //quals += read.getBaseQualityString().charAt(offset);
                //offsetString += i;
                //System.out.printf("  [%2d] [%s] %s%n", offset, read.getReadString().charAt(offset), read.getReadString());
        } 

        GenotypeEvidence all = new GenotypeEvidence(bases, ref);

        // P(q|G) - prob of nonref mixture given the genotype
        float qobs = all.q; // observed percent of non-ref bases
        double G; // % non-ref bases in observed
        if (qobs >= 0.1) {
            all.print();
            System.out.format("q %.2f | ", all.q);
            System.out.format("%s | ", context.getLocation());
            System.out.format("Total %4d | ", context.numReads());
            System.out.println();
            for (int q = 0; q < all.allbases; q ++) {
                for (G = 0.01; G <= 1.0; G += 0.49) { // iterate over: ref (0%), het (50%) and hom (100%) nonref bases observed
                    //double pqG = binomialProb(all.allbases - all.refbases, all.allbases, G);
                    double pqG = binomialProb(q, all.allbases, G);
                        //all.print();
                    System.out.format("P(q|G) %.3f | ", pqG);
                }
                System.out.println();
            }
            long stop_tm = currentTimeMillis();
            System.out.format("%.3fs\n", (float)(stop_tm - start_tm) / 1000);
        }
        return 1;
    }

    static double binomialProb(int k, int n, double p) {
        // k - numebr of successes
        // n - number of Bernoulli trials
        // p - probability of success

        return (double)nchoosek(n, k) * Math.pow(p, k) * Math.pow(1-p, n-k);
    }

    static int nchoosek(int n, int k) {
        int t = 1;

        int m = n - k;
        if (k < m) {
            k = m;
        }

        for (int i = n, j = 1; i > k; i--, j++) {
            t = t * i / j;
        }

        return t;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
