package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMTag;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.Reads;
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
    protected static Logger logger = Logger.getLogger(ReadFormattingIterator.class);    

    /**
     * Iterator to which to pass
     */
    private StingSAMIterator wrappedIterator;

    /**
     * Decorate the given iterator inside a ReadWrappingIterator.
     * @param wrappedIterator iterator
     */
    public ReadFormattingIterator(StingSAMIterator wrappedIterator) {
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
        SAMRecord read = wrappedIterator.next();

        // if we don't have a read group, set one.
        // TODO: Straw poll to see whether this is really required.        
        if (read.getAttribute(SAMTag.RG.toString()) == null && read.getReader() != null) {
            List<SAMReadGroupRecord> readGroups = read.getReader().getFileHeader().getReadGroups();
            if (readGroups.size() == 1) {
                read.setAttribute(SAMTag.RG.toString(), readGroups.get(0).getReadGroupId());
                read.setAttribute(SAMTag.SM.toString(), readGroups.get(0).getReadGroupId());
            } else {
                logger.warn("Unable to set read group of ungrouped read: unable to pick default group, there are " + readGroups.size() + " possible.");
            }
        }

        return new GATKSAMRecord(read);
    }

    /**
     * Remove the current element from the list.  Unsupported in this wrapper.
     */
    public void remove() { throw new UnsupportedOperationException("Cannot remove from a ReadWrappingIterator"); }
}
