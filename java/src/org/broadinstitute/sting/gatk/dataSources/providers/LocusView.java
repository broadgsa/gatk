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
import java.util.NoSuchElementException;
import java.util.Collection;
import java.util.Collections;
import java.util.Arrays;

import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;
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

public abstract class LocusView extends LocusContextIterator implements View {
    /**
     * The shard bounding this view.
     */
    protected Shard shard;

    /**
     * Source info for this view.  Informs the class about downsampling requirements.
     */
    private Reads sourceInfo;

    /**
     * The actual locus context iterator.
     */
    private LocusContextIterator loci;

    /**
     * The next locus context from the iterator.  This value must always be within
     * the shard; if its null, there's nothing for the consumer to look at. 
     */
    private LocusContext nextLocusContext = null;

    public LocusView(ShardDataProvider provider) {
        this.shard = provider.getShard();
        
        Iterator<SAMRecord> reads = new FilteringIterator(provider.getReadIterator(), new LocusStreamFilterFunc());
        this.sourceInfo = provider.getReadIterator().getSourceInfo();

        this.loci = new LocusContextIteratorByHanger(reads);
        seedNextLocusContext();

        provider.register(this);
    }

    /**
     * Only one view of the locus is supported at any given time.
     * @return A list consisting of all other locus views.
     */
    public Collection<Class<? extends View>> getConflictingViews() {
        return Arrays.<Class<? extends View>>asList(LocusView.class,ReadView.class);
    }

    /**
     * Close this view.
     */
    public void close() {
        // Set everything to null with the hope of failing fast.
        shard = null;
        sourceInfo = null;
        loci = null;

        super.close();
    }

    /**
     * Is there another covered locus context bounded by this view.
     * @return True if another covered locus context exists.  False otherwise.
     */
    public abstract boolean hasNext();

    /**
     * Returns the next covered locus context in the shard.
     * @return Next covered locus context in the shard.
     * @throw NoSuchElementException if no such element exists.
     */
    public abstract LocusContext next();

    /**
     * Unsupported.
     * @throw UnsupportedOperationException always.
     */
    public void remove() {
        throw new UnsupportedOperationException("Unable to remove elements from this queue.");
    }

    /**
     * Is there another locus context bounded by this shard.
     * @return True if another locus context is bounded by this shard.
     */
    protected boolean hasNextLocusContext() {
        return nextLocusContext != null && !nextLocusContext.getLocation().isPast(shard.getGenomeLoc());
    }

    /**
     * Get the next locus context bounded by this shard.
     * @return Next locus context bounded by this shard.
     * @throw NoSuchElementException if the next element is missing.
     */
    protected LocusContext nextLocusContext() {
        if( nextLocusContext == null || nextLocusContext.getLocation().isPast(shard.getGenomeLoc()) )
            throw new NoSuchElementException("No more elements remain in locus context queue.");

        // Cache the current and apply filtering.
        LocusContext current = nextLocusContext;

        // Find the next.
        if( loci.hasNext() ) {
            nextLocusContext = loci.next();
            if( sourceInfo.getDownsampleToCoverage() != null )
                current.downsampleToCoverage( sourceInfo.getDownsampleToCoverage() );                                 
            if( nextLocusContext.getLocation().isPast(shard.getGenomeLoc()) )
                nextLocusContext = null;
        }
        else
            nextLocusContext = null;

        return current;
    }

    /**
     * Seed the nextLocusContext variable with the contents of the next locus context (if one exists).
     */
    private void seedNextLocusContext() {
        if( loci.hasNext() )
            nextLocusContext = loci.next();

        // Iterate past cruft at the beginning to the first locus in the shard.
        while( nextLocusContext != null && nextLocusContext.getLocation().isBefore(shard.getGenomeLoc()) && loci.hasNext() )
            nextLocusContext = loci.next();

        // If nothing in the shard was found, indicate that by setting nextLocusContext to null.
        if( nextLocusContext != null && nextLocusContext.getLocation().isBefore(shard.getGenomeLoc()) )
            nextLocusContext = null;
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
