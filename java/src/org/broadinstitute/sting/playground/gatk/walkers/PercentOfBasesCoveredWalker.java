package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.cmdLine.Argument;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: Nov 25, 2009
 * Time: 7:54:42 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
public class PercentOfBasesCoveredWalker extends LocusWalker<Boolean, Pair<Integer,Integer>> implements TreeReducible<Pair<Integer,Integer>> {
    @Argument(fullName="minimumDepth", shortName="minDepth", doc="Minimum depth beyond which to discount the base as uncovered", required=false)
    int minDepth = 0;
    @Argument(fullName="includeName",doc="include this name in the output of this walker (e.g. NAME: percent of bases covered . . .) ", required=false)
    String name = null;
    public void initialize() {
        // do nothing
    }

    public Pair<Integer,Integer> reduceInit() {
        return new Pair<Integer,Integer>(0,1); // robust initialization -- seen one base with zero coverage
    }

    public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return context.size() > minDepth;
    }

    public Pair<Integer,Integer> reduce(Boolean map, Pair<Integer,Integer> prevReduce) {
        if ( map.booleanValue() ) {
            return new Pair<Integer,Integer>(prevReduce.first+1,prevReduce.second+1);
        } else {
            return new Pair<Integer,Integer>(prevReduce.first,prevReduce.second+1);
        }
    }

    public Pair<Integer,Integer> treeReduce( Pair<Integer,Integer> left, Pair<Integer,Integer> right) {
        right.set(right.first+left.first,right.second+left.second);
        return right;
    }

    public void onTraversalDone(Pair<Integer,Integer> reduceFinal) {
        String msg;
        if (name == null)
            msg = "Percent of bases covered at "+minDepth+"x=";
        else
            msg = name+": Percent of bases covered at "+minDepth+"x=";
        
        out.printf("%s%s",msg,pair2percent(reduceFinal));
    }

    private String pair2percent(Pair<Integer,Integer> p ) {
        return String.format("%d%s", (int) Math.floor( ( ( double) p.first*100.0) / p.second ), "%" );
    }
}
