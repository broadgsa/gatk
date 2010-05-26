package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Arrays;
import java.util.Collection;
import java.util.NoSuchElementException;

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

public abstract class LocusView extends LocusIterator implements View {
    /**
     * The locus bounding this view.
     */
    protected GenomeLoc locus;

    /**
     * Source info for this view.  Informs the class about downsampling requirements.
     */
    private Reads sourceInfo;

    /**
     * The actual locus context iterator.
     */
    private LocusIterator loci;

    /**
     * The next locus context from the iterator.  This value must always be within
     * the shard; if its null, there's nothing for the consumer to look at. 
     */
    private AlignmentContext nextLocus = null;

    public LocusView(LocusShardDataProvider provider) {
        this.locus = provider.getLocus();
        
        this.sourceInfo = provider.getSourceInfo();
        this.loci = provider.getLocusIterator();

        seedNextLocus();

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
        locus = null;
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
    public abstract AlignmentContext next();

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
    protected boolean hasNextLocus() {
        return nextLocus != null;
    }

    /**
     * Get the next locus context bounded by this shard.
     * @return Next locus context bounded by this shard.
     * @throw NoSuchElementException if the next element is missing.
     */
    protected AlignmentContext nextLocus() {
        if(nextLocus == null)
            throw new NoSuchElementException("No more elements remain in locus context queue.");

        // Cache the current and apply filtering.
        AlignmentContext current = nextLocus;

        // Find the next.
        seedNextLocus();
        if( sourceInfo.getDownsamplingMethod().type == DownsampleType.ALL_READS && sourceInfo.getDownsamplingMethod().toCoverage != null )
            current.downsampleToCoverage( sourceInfo.getDownsamplingMethod().toCoverage );
        
        // if the current loci isn't null, get the overflow tracker and pass it to the alignment context
        if ((this.loci != null))
            current.setLocusOverflowTracker(loci.getLocusOverflowTracker());
        return current;
    }

    /**
     * Seed the nextLocus variable with the contents of the next locus (if one exists).
     */
    private void seedNextLocus() {
        //System.out.printf("loci is %s%n", loci);
        if( !loci.hasNext() ) {
            nextLocus = null;
            return;
        }

        nextLocus = loci.next();

        // If the location of this shard is available, trim the data stream to match the shard.
        // TODO: Much of this functionality is being replaced by the WindowMaker.
        if(locus != null) {
            // Iterate through any elements not contained within this shard.
            while( nextLocus != null && !isContainedInShard(nextLocus.getLocation()) && loci.hasNext() )
                nextLocus = loci.next();

            // If nothing in the shard was found, indicate that by setting nextAlignmentContext to null.
            if( nextLocus != null && !isContainedInShard(nextLocus.getLocation()) )
                nextLocus = null;
        }
    }

    /**
     * Is this location contained in the given shard.
     * @param location Location to check.
     * @return True if the given location is contained within the shard.  False otherwise.
     */
    private boolean isContainedInShard(GenomeLoc location) {
        return locus.containsP(location);
    }
}
