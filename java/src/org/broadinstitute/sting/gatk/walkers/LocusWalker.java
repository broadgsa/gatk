package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.READS)
@Requires({DataSource.READS,DataSource.REFERENCE, DataSource.REFERENCE_BASES})
public abstract class LocusWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    public abstract MapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context);
}
