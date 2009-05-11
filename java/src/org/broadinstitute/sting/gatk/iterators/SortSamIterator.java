package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.ComparableSAMRecord;

import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

public class SortSamIterator implements Iterator<SAMRecord> {

    Iterator<ComparableSAMRecord> it;

    public SortSamIterator(Iterator<SAMRecord> unsortedIter, int maxSorts) {

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

    public boolean hasNext() { return it.hasNext(); }
    public SAMRecord next()  { return it.next().getRecord(); }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }
}
