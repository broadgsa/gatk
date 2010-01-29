package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.*;

import java.util.TreeSet;
import java.util.Iterator;
import java.util.List;

/**
 * @author ebanks
 * SortingSAMFileWriter
 *
 * this class extends the samtools SAMFileWriter class and caches reads for N loci so that reads
 * can be emitted out of order (provided they are within the N-locus window)
 *
 */
public class SortingSAMFileWriter implements SAMFileWriter {

    // the base writer from Picard
    private SAMFileWriter baseWriter;

    // the window over which we agree to accumulate reads
    private int window;

    // the reads we are accumulating
    private TreeSet<SAMRecord> cachedReads = new TreeSet<SAMRecord>(new SAMRecordCoordinateComparatorWithUnmappedReads());


    /**
     * Constructor
     *
     * @param baseWriter   the real SAMFileWriter
     * @param window       the window over which we agree to store reads
     */
    public SortingSAMFileWriter(SAMFileWriter baseWriter, int window) {
        this.baseWriter = baseWriter;
        this.window = window;
    }

    /**
     * Add a read to the writer for emission
     *
     * @param read   the read to emit
     */
    public void addAlignment(SAMRecord read) {
        // at a new contig, clear the cache
        if ( cachedReads.size() > 0 && cachedReads.first().getReferenceIndex() < read.getReferenceIndex() )
            clearCache();

        long currentPos = read.getAlignmentStart();

        Iterator<SAMRecord> iter = cachedReads.iterator();
        while ( iter.hasNext() ) {
            SAMRecord cachedRead = iter.next();
            if ( currentPos - cachedRead.getAlignmentStart() >= window ) {
                baseWriter.addAlignment(cachedRead);
                iter.remove();                
            } else {
                break;
            }
        }

        cachedReads.add(read);
    }

    /**
     * Add a list of reads to the writer for emission; the reads do NOT need to be sorted
     *
     * @param reads   the reads to emit
     */
    public void addAlignments(List<SAMRecord> reads) {
        if ( reads.size() == 0 )
            return;

        // at a new contig, clear the cache
        if ( cachedReads.size() > 0 && cachedReads.first().getReferenceIndex() < reads.get(0).getReferenceIndex() )
            clearCache();

        cachedReads.addAll(reads);

        // get the last read in the cache
        SAMRecord last = cachedReads.last();

        long currentPos = last.getAlignmentStart();

        Iterator<SAMRecord> iter = cachedReads.iterator();
        while ( iter.hasNext() ) {
            SAMRecord cachedRead = iter.next();
            if ( currentPos - cachedRead.getAlignmentStart() >= window ) {
                baseWriter.addAlignment(cachedRead);
                iter.remove();
            } else {
                break;
            }
        }
    }

    /**
     * get the SAM file header
     */
    public SAMFileHeader getFileHeader() {
        return baseWriter.getFileHeader();
    }

    /**
     * close this writer by clearing the cache
     */
    public void close() {
        clearCache();
    }

    private void clearCache() {
        Iterator<SAMRecord> iter = cachedReads.iterator();
        while ( iter.hasNext() )
            baseWriter.addAlignment(iter.next());
        cachedReads.clear();
    }
}
