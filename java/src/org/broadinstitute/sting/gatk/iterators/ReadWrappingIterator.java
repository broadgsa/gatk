package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * An iterator which wraps each SAMRecord inside a wrapper class, bringing new functionality to the read while
 * presenting the original SAMRecord interface.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadWrappingIterator implements StingSAMIterator {
    /**
     * Iterator to which to pass
     */
    private StingSAMIterator wrappedIterator;

    /**
     * Decorate the given iterator inside a ReadWrappingIterator.
     * @param wrappedIterator iterator
     */
    public ReadWrappingIterator(StingSAMIterator wrappedIterator) {
        this.wrappedIterator = wrappedIterator;
    }

    /**
     * Get metadata about the reads' sources, etc.
     * @return Source info about the reads.
     */
    public Reads getSourceInfo() {
        return wrappedIterator.getSourceInfo();
    }

    /**
     * Convenience function for use in foreach loops.  Dangerous because it does not actually
     * reset the iterator.
     * @return An iterator through the current data stream.
     */
    public StingSAMIterator iterator() {
        // NOTE: this iterator doesn't perform any kind of reset operation; it just returns itself.
        //       can we do something better?  Do we really have to provide support for the Iterable interface?
        return this;
    }

    /**
     * Close this iterator.
     */
    public void close() {
        wrappedIterator.close();
    }

    /**
     * Does the iterator contain more values?
     * @return True if there are more left to return, false otherwise.
     */
    public boolean hasNext() {
        return wrappedIterator.hasNext();
    }

    /**
     * Get the next value in the sequence.
     * @return Next value in the sequence.  By convention, a NoSuchElementException should be thrown if
     *         no next exists.
     */
    public SAMRecord next() {
        return new GATKSAMRecord(wrappedIterator.next());
    }

    /**
     * Remove the current element from the list.  Unsupported in this wrapper.
     */
    public void remove() { throw new UnsupportedOperationException("Cannot remove from a ReadWrappingIterator"); }
}
