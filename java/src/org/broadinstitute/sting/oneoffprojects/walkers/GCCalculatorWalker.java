package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.commandline.Output;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.io.PrintStream;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date May 19, 2010
 */
public class GCCalculatorWalker extends RefWalker<Pair<Set<GenomeLoc>,Boolean>, Map<GenomeLoc,Pair<Long,Long>>> {
    @Output
    PrintStream out;

    public Map<GenomeLoc,Pair<Long,Long>> reduceInit() {
        return new HashMap<GenomeLoc,Pair<Long,Long>>();
    }

    public Pair<Set<GenomeLoc>,Boolean> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || tracker.getReferenceMetaData("interval_list") == null ) {
            return null;
        } else {
            Set<GenomeLoc> overlappingIntervals = new HashSet<GenomeLoc>();
            for ( GATKFeature f : tracker.getGATKFeatureMetaData("interval_list",true) ) {
                overlappingIntervals.add( f.getLocation() );
            }

            return new Pair<Set<GenomeLoc>,Boolean>(overlappingIntervals, ref.getBaseIndex() == BaseUtils.cIndex || ref.getBaseIndex() == BaseUtils.gIndex );
        }
    }

    public Map<GenomeLoc,Pair<Long,Long>> reduce(Pair<Set<GenomeLoc>,Boolean> map, Map<GenomeLoc,Pair<Long,Long>> prevReduce) {
        if ( map == null ) {
            return prevReduce;
        }

        for ( GenomeLoc loc : map.first ) {
            if ( ! prevReduce.keySet().contains(loc) ) {
                prevReduce.put(loc,new Pair<Long,Long>(0l,0l));
            }

            prevReduce.get(loc).first ++;
            if ( map.second ) {
                prevReduce.get(loc).second ++;
            }
        }

        return prevReduce;
    }

    public void onTraversalDone(Map<GenomeLoc,Pair<Long,Long>> reduced ) {
        for ( Map.Entry<GenomeLoc,Pair<Long,Long>> gcCounts : reduced.entrySet() ) {
            double gc_content = ( (double) gcCounts.getValue().second )/( (double) gcCounts.getValue().first );
            out.printf("%s\t%.2f%n",gcCounts.getKey().toString(),100*gc_content);
        }
    }
}
