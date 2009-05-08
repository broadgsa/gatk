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
public abstract class ReadWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    public boolean requiresOrderedReads() { return false; }
    
    // Do we actually want to operate on the context?
    /** Must return true for reads that need to be processed. Reads, for which this method return false will
     * be skipped by the engine and never passed to the walker.
     */
    public boolean filter(LocusContext context, SAMRecord read) {
        // We are keeping all the reads
        return true;
    }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public abstract MapType map(LocusContext context, SAMRecord read);
}
