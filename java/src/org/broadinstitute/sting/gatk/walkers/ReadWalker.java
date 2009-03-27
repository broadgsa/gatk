package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReadWalker<MapType, ReduceType> extends Walker<ReduceType> {
    public boolean requiresOrderedReads() { return false; }
    
    // Do we actually want to operate on the context?
    public boolean filter(LocusContext context, SAMRecord read) {
        // We are keeping all the reads
        return true;
    }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public abstract MapType map(LocusContext context, SAMRecord read);

    // Given result of map function
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);
}
