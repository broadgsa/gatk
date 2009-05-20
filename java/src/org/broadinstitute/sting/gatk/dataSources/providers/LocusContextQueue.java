package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.LocusContextIteratorByHanger;
import org.broadinstitute.sting.gatk.iterators.LocusContextIterator;
import org.broadinstitute.sting.gatk.traversals.TraversalStatistics;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.samtools.SAMRecord;

import java.util.Iterator;

import edu.mit.broad.picard.filter.FilteringIterator;
import edu.mit.broad.picard.filter.SamRecordFilter;
/**
 * User: hanna
 * Date: May 13, 2009
 * Time: 3:30:16 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A queue of locus context entries.
 */

public abstract class LocusContextQueue {
    protected Shard shard;

    private Reads sourceInfo;
    private LocusContextIterator loci;    

    public LocusContextQueue(ShardDataProvider provider) {
        Iterator<SAMRecord> reads = new FilteringIterator(provider.getReadIterator(), new LocusStreamFilterFunc());
        this.loci = new LocusContextIteratorByHanger(reads);
        this.sourceInfo = provider.getReadIterator().getSourceInfo();
        this.shard = provider.getShard();
    }


    /**
     * Get the locus context at the given position.
     * @return Locus context, or null if no locus context exists at this position.
     */
    public abstract LocusContext peek();

    /**
     * Seek to the given point the queue of locus contexts.
     * @param target Target base pair to which to seek.  Must be a single base pair.
     * @return an instance of itself for parameter chaining.
     */
    public abstract LocusContextQueue seek(GenomeLoc target);

    /**
     * Gets the next locus context, applying filtering as necessary.
     * @return Locus context to work with.
     */
    protected LocusContext getNextLocusContext() {
        LocusContext next = loci.next();
        if( sourceInfo.getDownsampleToCoverage() != null )
            next.downsampleToCoverage( sourceInfo.getDownsampleToCoverage() );
        return next;
    }

    /**
     * hasNext()-style iterator for base iterator. 
     * @return
     */
    protected boolean hasNextLocusContext() {
        return loci.hasNext();
    }

    /**
     * Class to filter out un-handle-able reads from the stream.  We currently are skipping
     * unmapped reads, non-primary reads, unaligned reads, and duplicate reads.
     */
    private static class LocusStreamFilterFunc implements SamRecordFilter {
        SAMRecord lastRead = null;
        public boolean filterOut(SAMRecord rec) {
            boolean result = false;
            String why = "";
            if (rec.getReadUnmappedFlag()) {
                TraversalStatistics.nUnmappedReads++;
                result = true;
                why = "Unmapped";
            } else if (rec.getNotPrimaryAlignmentFlag()) {
                TraversalStatistics.nNotPrimary++;
                result = true;
                why = "Not Primary";
            } else if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                TraversalStatistics.nBadAlignments++;
                result = true;
                why = "No alignment start";
            } else if (rec.getDuplicateReadFlag()) {
                TraversalStatistics.nDuplicates++;
                result = true;
                why = "Duplicate reads";
            }
            else {
                result = false;
            }

            if (result) {
                TraversalStatistics.nSkippedReads++;
                //System.out.printf("  [filter] %s => %b %s", rec.getReadName(), result, why);
            } else {
                TraversalStatistics.nReads++;
            }
            return result;
        }
    }
}
