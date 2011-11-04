package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires({DataSource.READS, DataSource.REFERENCE_BASES})
@PartitionBy(PartitionType.READ)
public abstract class ReadWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    public boolean requiresOrderedReads() { return false; }
    
    // Do we actually want to operate on the context?
    /** Must return true for reads that need to be processed. Reads, for which this method return false will
     * be skipped by the engine and never passed to the walker.
     */
    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        // We are keeping all the reads
        return true;
    }

    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    public abstract MapType map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker);
}
