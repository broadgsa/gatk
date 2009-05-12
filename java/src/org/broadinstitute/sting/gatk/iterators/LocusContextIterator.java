package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.util.CloseableIterator;

import java.util.Iterator;

import org.broadinstitute.sting.gatk.LocusContext;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 */
public abstract class LocusContextIterator implements Iterable<LocusContext>, CloseableIterator<LocusContext> {
    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public Iterator<LocusContext> iterator() {
        return this;
    }

    public void close() {
        //this.it.close();
    }

    public abstract boolean hasNext();
    public abstract LocusContext next();

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }
}
