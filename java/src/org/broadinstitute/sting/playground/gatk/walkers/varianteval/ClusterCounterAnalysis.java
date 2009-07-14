package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.IntervalRod;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashSet;
import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Jun 4, 2009
 * Time: 4:38:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class ClusterCounterAnalysis extends BasicVariantAnalysis {
    ArrayList<HashSet<GenomeLoc>> variantsWithClusters;
    int[] neighborWiseBoundries = {1, 2, 5, 10, 20, 50, 100};
    AllelicVariant lastVariant = null;
    GenomeLoc lastVariantInterval = null;
    int nSeen = 0;

    public ClusterCounterAnalysis() {
        super("cluster_counter_analysis");

        variantsWithClusters = new ArrayList<HashSet<GenomeLoc>>(neighborWiseBoundries.length);
        for (int i = 0; i < neighborWiseBoundries.length; i++)
            variantsWithClusters.add(new HashSet<GenomeLoc>());
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
                        nSeen++;
                        StringBuilder s = new StringBuilder();
                        for (int i = 0; i < neighborWiseBoundries.length; i++) {
                            int maxDist = neighborWiseBoundries[i];
                            s.append(String.format("%d ", d <= maxDist ? maxDist : 0));
                            if ( d <= maxDist ) {
                                variantsWithClusters.get(i).add(eL);
                            }
                        }
                        r = String.format("snp_within_cluster %d %s %s %s", d, eL, lvL, s.toString());
                    }
                }
            }

            lastVariant = eval;
            lastVariantInterval = interval;
        }

        return r;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("snps_counted_for_neighbor_distances %d", nSeen));
        s.add(String.format("description        maxDist count"));
        for ( int i = 0; i < neighborWiseBoundries.length; i++ ) {
            int maxDist = neighborWiseBoundries[i];
            int count = variantsWithClusters.get(i).size();
            s.add(String.format("snps_within_clusters_of_size %10d %10d", maxDist, count));
        }

        return s;
    }
}