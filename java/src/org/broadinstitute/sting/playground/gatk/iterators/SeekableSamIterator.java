package org.broadinstitute.sting.playground.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;

import java.util.Iterator;

import org.broadinstitute.sting.playground.gatk.iterators.SeekableSamIteration;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:24:38 AM
 * To change this template use File | Settings | File Templates.
 */
public class SeekableSamIterator implements Iterator<SAMRecord>, SeekableSamIteration {
    protected Iterator<SAMRecord> it;
    protected SAMFileReader reader;

    public SeekableSamIterator(Iterator<SAMRecord> it, SAMFileReader reader) {
        this.it = it;
        this.reader = reader;
    }

    public boolean supportsSeeking() { return true; }

    public void queryOverlapping( final String contig, final int start, final int stop ) {
        this.it = reader.queryOverlapping( contig, start, stop );
    }

    public void query(final String contig, final int start, final int stop, final boolean contained) {
        this.it = reader.query( contig, start, stop, contained );
    }

    public void queryContained(final String contig, final int start, final int stop) {
        this.it = reader.queryContained( contig, start, stop );
    }

    public boolean hasNext() { return it.hasNext(); }
    public SAMRecord next() { return it.next(); }
    public void remove () { it.remove(); }
}