package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMTag;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.apache.log4j.Logger;

import java.util.List;

/**
 * An iterator which does post-processing of a read, including potentially wrapping
 * the read in something with a compatible interface or replacing the read entirely.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadFormattingIterator implements StingSAMIterator {
    /**
     * Logger.
     */
    final protected static Logger logger = Logger.getLogger(ReadFormattingIterator.class);    

    /**
     * Iterator to which to pass
     */
    private StingSAMIterator wrappedIterator;

    /**
     * True if original base qualities should be used.
     */
    private final boolean useOriginalBaseQualities;

    /**
     * Decorate the given iterator inside a ReadWrappingIterator.
     * @param wrappedIterator iterator
     * @param useOriginalBaseQualities true if original base qualities should be used
     */
    public ReadFormattingIterator(StingSAMIterator wrappedIterator, boolean useOriginalBaseQualities) {
        this.wrappedIterator = wrappedIterator;
        this.useOriginalBaseQualities = useOriginalBaseQualities;
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
        return new GATKSAMRecord(wrappedIterator.next(), useOriginalBaseQualities);
    }

    /**
     * Remove the current element from the list.  Unsupported in this wrapper.
     */
    public void remove() { throw new UnsupportedOperationException("Cannot remove from a ReadWrappingIterator"); }
}
