package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.BasicPileup;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 28, 2009
 * Time: 3:19:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class IndelCounterWalker extends LocusWalker<Integer,Integer> {

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce( Integer prevReduce, Integer map ) {
        return map + prevReduce;
    }

    public Integer treeReduce( Integer reduce1, Integer reduce2 ) {
        return reduce(reduce1,reduce2);
    }

    public Integer map ( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        if (BaseUtils.isRegularBase(ref.getBase()) ) {
            return numIndels(context);
        } else {
            return 0;
        }
    }

    public Integer numIndels ( AlignmentContext context ) {
        String[] indelPileup = BasicPileup.indelPileup(context.getReads(),context.getOffsets());
        // number of indels is the number of non-"null" indeces
        int nIndel = 0;
        for( String indel : indelPileup ) {
            if ( ! indel.equals("null") ) {  // currently how non-indel bases are represented in pileup
                nIndel ++;
            }
        }

        return nIndel;
    }
}
