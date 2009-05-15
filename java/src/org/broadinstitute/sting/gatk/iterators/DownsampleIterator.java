package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;

import java.util.Iterator;
import java.util.Random;

import org.broadinstitute.sting.gatk.Reads;


public class DownsampleIterator implements StingSAMIterator {

    StingSAMIterator it;
    Random generator;
    int cutoff;
    SAMRecord next;

    public DownsampleIterator(StingSAMIterator it, double fraction) {
        this.it = it;
        generator = new Random();
        cutoff = (int)(fraction * 10000);
        next = getNextRecord();
    }

    /**
     * Retrieves information about reads sources.
     * @return Info about the sources of reads.
     */
    public Reads getSourceInfo() {
        return it.getSourceInfo();
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
            if ( generator.nextInt(10000) < cutoff )
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