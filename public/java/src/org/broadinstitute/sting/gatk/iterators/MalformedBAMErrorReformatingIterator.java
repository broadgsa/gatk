package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * Traps BAM formatting errors in underlying iterator and rethrows meaningful GATK UserExceptions
 */
public class MalformedBAMErrorReformatingIterator implements CloseableIterator<SAMRecord> {
    File source;
    CloseableIterator<SAMRecord> it;

    public MalformedBAMErrorReformatingIterator(final File source, final CloseableIterator<SAMRecord> it) {
        this.it = it;
        this.source = source;
    }

    public boolean hasNext() {
        try {
            return this.it.hasNext();
        } catch ( RuntimeException e ) { // we need to catch RuntimeExceptions here because the Picard code is throwing them (among SAMFileExceptions) sometimes
            throw new UserException.MalformedBAM(source, e.getMessage());
        }
    }

    public SAMRecord next() {
        try {
            return it.next();
        } catch ( RuntimeException e ) { // we need to catch RuntimeExceptions here because the Picard code is throwing them (among SAMFileExceptions) sometimes
            throw new UserException.MalformedBAM(source, e.getMessage());
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    public void close() { it.close(); }
    public Iterator<SAMRecord> iterator() { return this; }
}
