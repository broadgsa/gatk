package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Apr 23, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires({DataSource.READS,DataSource.REFERENCE})
public abstract class LocusWindowWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, String ref, LocusContext context) {
        return true;    // We are keeping all the intervals
    }

    // do we care about reads that are not part of our intervals?
    public boolean actOnNonIntervalReads() {
        return false;   // Don't act on them
    }

    // What do we do with the reads that are not part of our intervals?
    public void nonIntervalReadAction(SAMRecord read) { }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public abstract MapType map(RefMetaDataTracker tracker, String ref, LocusContext context);

    // Given result of map function
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);
}