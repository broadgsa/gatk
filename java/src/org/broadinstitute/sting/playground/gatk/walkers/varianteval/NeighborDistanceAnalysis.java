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
import java.io.File;
import java.io.FileNotFoundException;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Jun 4, 2009
 * Time: 4:38:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class NeighborDistanceAnalysis extends ViolationVariantAnalysis {
    ArrayList<Long> neighborWiseDistances;
    int[] neighborWiseBoundries = {1, 2, 5, 10, 20, 50, 100, 1000, 10000};

    AllelicVariant lastVariant = null;
    GenomeLoc lastVariantInterval = null;
    PrintStream violationsOut = null;

    public NeighborDistanceAnalysis() {
        super("neighbor_distances");
        neighborWiseDistances = new ArrayList<Long>();
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
                        neighborWiseDistances.add(d);
                        r = String.format("neighbor-distance %d %s %s", d, eL, lvL);
                    }
                }
            }
            
            lastVariant = eval;
            lastVariantInterval = interval;
        }
        
        return r;
    }

    public List<String> done() {
        int[] pairCounts = new int[neighborWiseBoundries.length];
        Arrays.fill(pairCounts, 0);
        for ( long dist : neighborWiseDistances ) {
            boolean done = false;
            for ( int i = 0; i < neighborWiseBoundries.length && ! done ; i++ ) {
                int maxDist = neighborWiseBoundries[i];
                if ( dist <= maxDist ) {
                    pairCounts[i]++;
                    done = true;
                }
            }
        }

        List<String> s = new ArrayList<String>();
        s.add(String.format("snps_counted_for_neighbor_distances %d", neighborWiseDistances.size()));
        int cumulative = 0;
        s.add(String.format("description        maxDist count cumulative"));
        for ( int i = 0; i < neighborWiseBoundries.length; i++ ) {
            int maxDist = neighborWiseBoundries[i];
            int count = pairCounts[i];
            cumulative += count;
            s.add(String.format("snps_immediate_neighbors_within_bp %d  %d  %d", maxDist, count, cumulative));
        }

        return s;
    }
}
