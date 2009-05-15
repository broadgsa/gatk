package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.ComparableSAMRecord;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.Reads;

import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

// TODO: Deprecate? 
// I don't think we need this if we're only allowing sorted and indexed BAM Files in the GATK - Aaron
public class SortSamIterator implements StingSAMIterator {
    private Reads sourceInfo;
    private Iterator<ComparableSAMRecord> it;

    public SortSamIterator(StingSAMIterator unsortedIter, int maxSorts) {
        sourceInfo = unsortedIter.getSourceInfo();
        ArrayList<ComparableSAMRecord> list = new ArrayList<ComparableSAMRecord>();
        while (unsortedIter.hasNext()) {
            list.add(new ComparableSAMRecord(unsortedIter.next()));
            // limit how much can be sorted for now
            if (list.size() > maxSorts)
                throw new UnsupportedOperationException("Can not sort files with more than " + maxSorts + " reads on the fly!");
        }
        Collections.sort(list);
        it = list.iterator();
    }

    /**
     * Retrieves information about reads sources.
     * @return Info about the sources of reads.
     */
    public Reads getSourceInfo() {
        if( sourceInfo == null )
            throw new StingException("Unable to provide source info for the reads.  Please upgrade to the new data sharding framework.");
        return sourceInfo;
    }

    public boolean hasNext() { return it.hasNext(); }
    public SAMRecord next()  { return it.next().getRecord(); }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    public void close() {
        // nothing to do right now
    }

    public Iterator<SAMRecord> iterator() {
        return this;  
    }
}
