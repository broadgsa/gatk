package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.Iterator;

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
}
