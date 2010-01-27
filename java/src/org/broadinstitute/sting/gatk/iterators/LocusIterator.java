package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.util.CloseableIterator;

import java.util.Iterator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 */
public abstract class LocusIterator implements Iterable<AlignmentContext>, CloseableIterator<AlignmentContext> {
    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public Iterator<AlignmentContext> iterator() {
        return this;
    }

    public void close() {
        //this.it.close();
    }

    public abstract boolean hasNext();
    public abstract AlignmentContext next();

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    /**
     * a method for getting the overflow tracker, which is used to track sites at which the read count exceeds the
     * pile-up threshold set on the command line
     * 
     * @return the overflow tracker, null if no tracker exists
     */
    public abstract LocusOverflowTracker getLocusOverflowTracker();
}
