package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class DuplicateWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(GenomeLoc loc, byte[] refBases, LocusContext context, List<SAMRecord> duplicateReads) {
        return true;    // We are keeping all the reads
    }

    /**
     * These two functions state whether we're don't make any sense without reads (requiresRead())
     * or whether we can't take any reads at all (cannotHandleRead())
     */
    public boolean requiresReads()     { return true; }
    public boolean cannotHandleReads() { return false; }


    /**
     * Called by the traversal engine to decide whether to send non-duplicates as lists of
     * singleton reads to the map function.  By default it's false.
     * 
     * @return true if you want to see non duplicates during the traversal
     */
    public boolean mapUniqueReadsTooP() { return false; }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public abstract MapType map(GenomeLoc loc, byte[] refBases, LocusContext context, List<SAMRecord> duplicateReads);

    // Given result of map function
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);
}