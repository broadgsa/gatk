package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.util.Iterator;


public class DownsampleIterator implements StingSAMIterator {

    StingSAMIterator it;
    int cutoff;
    SAMRecord next;

    public DownsampleIterator(StingSAMIterator it, double fraction) {
        this.it = it;
        cutoff = (int)(fraction * 10000);
        next = getNextRecord();
    }

    public boolean hasNext() {
        return next != null;
    }

    public SAMRecord next()  {
        SAMRecord result = next;
        next = getNextRecord();
        return result;
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    private SAMRecord getNextRecord() {
        while ( true ) {
            if ( !it.hasNext() )
                return null;
            SAMRecord rec = it.next();
            if ( GenomeAnalysisEngine.getRandomGenerator().nextInt(10000) < cutoff )
                return rec;
        }
    }

    public void close() {
        it.close();
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}