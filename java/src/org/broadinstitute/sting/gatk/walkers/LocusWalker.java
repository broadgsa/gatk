package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.READS)
public abstract class LocusWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    /**
     * These two functions state whether we're don't make any sense without reads (requiresRead())
     * or whether we can't take any reads at all (cannotHandleRead())
     */
    public boolean requiresReads()     { return false; }
    public boolean cannotHandleReads() { return false; }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public abstract MapType map(RefMetaDataTracker tracker, char ref, LocusContext context);
}
