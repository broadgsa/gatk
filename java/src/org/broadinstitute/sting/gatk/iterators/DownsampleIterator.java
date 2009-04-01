package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;

import java.util.Random;
import java.util.Iterator;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;


public class DownsampleIterator implements Iterator<SAMRecord> {

    Iterator<SAMRecord> it;
    Random generator;
    int cutoff;
    SAMRecord next;

    public DownsampleIterator(Iterator<SAMRecord> it, double fraction) {
        this.it = it;
        generator = new Random();
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
            if ( generator.nextInt(10000) < cutoff )
                return rec;
        }
    }
}