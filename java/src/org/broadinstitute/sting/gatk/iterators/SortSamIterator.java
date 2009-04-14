package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 15, 2009
 * Time: 6:02:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class SortSamIterator implements Iterator<SAMRecord> {

    Iterator<ComparableSAMRecord> it;

    public SortSamIterator(Iterator<SAMRecord> unsortedIter, int maxSorts) {

        ArrayList<ComparableSAMRecord> list = new ArrayList<ComparableSAMRecord>();
        while (unsortedIter.hasNext()) {
            list.add(new ComparableSAMRecord(unsortedIter.next()));
            // choose an arbitrary length to limit sorting for now
            if (list.size() > maxSorts)
                throw new UnsupportedOperationException("Can not sort files with more than 100K reads on the fly!");
        }
        Collections.sort(list);
        it = list.iterator();
    }

    public boolean hasNext() { return it.hasNext(); }
    public SAMRecord next()  { return it.next().getRecord(); }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    private class ComparableSAMRecord implements Comparable<ComparableSAMRecord> {

        private SAMRecord record;

        public ComparableSAMRecord(SAMRecord record) {
            this.record = record;
        }

        public SAMRecord getRecord() {
            return record;
        }

        public int compareTo(ComparableSAMRecord o) {
            GenomeLoc myLoc = new GenomeLoc(record);
            GenomeLoc hisLoc = new GenomeLoc(o.getRecord());
            return myLoc.compareTo(hisLoc);
        }
    }
}
