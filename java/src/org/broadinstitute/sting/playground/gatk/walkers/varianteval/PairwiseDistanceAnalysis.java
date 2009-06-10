package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.IntervalRod;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Jun 4, 2009
 * Time: 4:38:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class PairwiseDistanceAnalysis extends BasicVariantAnalysis {
    ArrayList<Long> pairWiseDistances;
    int[] pairWiseBoundries = {1, 2, 5, 10, 20, 50, 100, 1000, 10000};

    AllelicVariant lastVariant = null;
    GenomeLoc lastVariantInterval = null;

    public PairwiseDistanceAnalysis() {
        super("pairwise_distances");
        pairWiseDistances = new ArrayList<Long>();
    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context) {
        String r = null;

        if ( eval != null ) {
            IntervalRod intervalROD = (IntervalRod)tracker.lookup("interval", null);
            GenomeLoc interval = intervalROD == null ? null : intervalROD.getLocation();

            if (lastVariant != null) {
                GenomeLoc eL = eval.getLocation();
                GenomeLoc lvL = lastVariant.getLocation();
                if (eL.getContigIndex() == lvL.getContigIndex()) {
                    long d = eL.distance(lvL);
                    if ( lastVariantInterval != null && lastVariantInterval.compareTo(interval) != 0) {
                        // we're on different intervals
                        //out.printf("# Excluding %d %s %s vs. %s %s%n", d, eL, interval, lvL, lastVariantInterval);
                    } else {
                        pairWiseDistances.add(d);
                        r = String.format("Pairwise-distance %d %s %s", d, eL, lvL);
                    }
                }
            }
            
            lastVariant = eval;
            lastVariantInterval = interval;
        }
        
        return r;
    }

    public List<String> done() {
        int[] pairCounts = new int[pairWiseBoundries.length];
        Arrays.fill(pairCounts, 0);
        for ( long dist : pairWiseDistances ) {
            boolean done = false;
            for ( int i = 0; i < pairWiseBoundries.length && ! done ; i++ ) {
                int maxDist = pairWiseBoundries[i];
                if ( dist <= maxDist ) {
                    pairCounts[i]++;
                    done = true;
                }
            }
        }

        List<String> s = new ArrayList<String>();
        s.add(String.format("snps counted for pairwise distance: %d", pairWiseDistances.size()));
        int cumulative = 0;
        for ( int i = 0; i < pairWiseBoundries.length; i++ ) {
            int maxDist = pairWiseBoundries[i];
            int count = pairCounts[i];
            cumulative += count;
            s.add(String.format("snps within %8d bp:    %d  %d", maxDist, count, cumulative));
        }

        return s;
    }
}
