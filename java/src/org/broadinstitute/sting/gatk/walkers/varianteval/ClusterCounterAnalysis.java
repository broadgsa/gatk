package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.IntervalRod;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.ArrayList;
import java.util.List;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
public class ClusterCounterAnalysis extends BasicVariantAnalysis implements GenotypeAnalysis, PopulationAnalysis {
    int minDistanceForFlagging = 5;

    /**
     * Threshold distances between neighboring SNPs we will track and assess
     */
    int[] neighborWiseBoundries = {1, 2, 5, 10, 20, 50, 100};
    int[] variantsWithClusters = new int[neighborWiseBoundries.length];
    /**
     * Keep track of the last variation we saw in the stream
     */
    Variation lastVariation = null;

    /**
     * The interval that the last variation occurred in.  Don't think that SNPs from different intervals are
     * too close together in hybrid selection.
     */
    GenomeLoc lastVariantInterval = null;
    int nSeen = 0;



    public ClusterCounterAnalysis() {
        super("cluster_counter_analysis");
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        String r = null;

        if ( eval != null && eval.isSNP() ) {
            IntervalRod intervalROD = tracker.lookup("interval",IntervalRod.class);
            GenomeLoc interval = intervalROD == null ? null : intervalROD.getLocation();

            if (lastVariation != null) {
                GenomeLoc eL = eval.getLocation();
                GenomeLoc lvL = lastVariation.getLocation();

                // if we are on the same contig, and we in the same interval
                if (eL.getContigIndex() == lvL.getContigIndex() && ! (lastVariantInterval != null && lastVariantInterval.compareTo(interval) != 0)) {
                    long d = eL.distance(lvL);
                    nSeen++;
                    StringBuilder s = new StringBuilder();
                    for (int i = 0; i < neighborWiseBoundries.length; i++) {
                        int maxDist = neighborWiseBoundries[i];
                        s.append(String.format("%d ", d <= maxDist ? maxDist : 0));
                        if ( d <= maxDist ) {
                            variantsWithClusters[i]++;
                        }
                    }

                    // lookup in master for performance reasons
                    if ( d <= minDistanceForFlagging && getMaster().includeViolations() )
                        r = String.format("snp_within_cluster %d %s %s %s", d, eL, lvL, s.toString());
                }
            }

            // eval is now the last variation
            lastVariation = eval;
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
            int count = variantsWithClusters[i];
            s.add(String.format("snps_within_clusters_of_size %10d %10d", maxDist, count));
        }

        return s;
    }
}