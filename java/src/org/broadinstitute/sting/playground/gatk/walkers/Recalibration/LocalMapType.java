package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 9, 2009
 * Time: 11:39:34 AM
 * To change this template use File | Settings | File Templates.
 */
class LocalMapType {

    public AlignmentContext context;
    public ReferenceContext ref;
    public RefMetaDataTracker tracker;

    public LocalMapType(AlignmentContext context, ReferenceContext ref, RefMetaDataTracker tracker) {
        this.context = context;
        this.ref = ref;
        this.tracker = tracker;
    }

    public int numReads() {
        return context.numReads();
    }

    public int qScore(int read) {
        return (int) context.getReads().get(read).getBaseQualities()[context.getOffsets().get(read)];
    }

}
