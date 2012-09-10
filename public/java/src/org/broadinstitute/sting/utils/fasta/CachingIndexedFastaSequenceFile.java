/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.fasta;

import net.sf.picard.PicardException;
import net.sf.picard.reference.FastaSequenceIndex;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Priority;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;

/**
 * A caching version of the IndexedFastaSequenceFile that avoids going to disk as often as the raw indexer.
 *
 * Thread-safe!  Uses a thread-local cache
 */
public class CachingIndexedFastaSequenceFile extends IndexedFastaSequenceFile {
    protected static final org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(CachingIndexedFastaSequenceFile.class);

    /** do we want to print debugging information about cache efficiency? */
    private static final boolean PRINT_EFFICIENCY = false;

    /** If we are printing efficiency info, what frequency should we do it at? */
    private static final int PRINT_FREQUENCY = 10000;

    /** The default cache size in bp */
    public static final long DEFAULT_CACHE_SIZE = 1000000;

    /** The cache size of this CachingIndexedFastaSequenceFile */
    final long cacheSize;

    /** When we have a cache miss at position X, we load sequence from X - cacheMissBackup */
    final long cacheMissBackup;

    // information about checking efficiency
    long cacheHits = 0;
    long cacheMisses = 0;

    /** Represents a specific cached sequence, with a specific start and stop, as well as the bases */
    private static class Cache {
        long start = -1, stop = -1;
        ReferenceSequence seq = null;
    }

    /**
     * Thread local cache to allow multi-threaded use of this class
     */
    private ThreadLocal<Cache> cache;
    {
        cache = new ThreadLocal<Cache> () {
            @Override protected Cache initialValue() {
                return new Cache();
            }
        };
    }

    /**
     * Same as general constructor but allows one to override the default cacheSize
     *
     * @param fasta
     * @param index
     * @param cacheSize
     */
    public CachingIndexedFastaSequenceFile(final File fasta, final FastaSequenceIndex index, final long cacheSize) {
        super(fasta, index);
        if ( cacheSize < 0 ) throw new IllegalArgumentException("cacheSize must be > 0");
        this.cacheSize = cacheSize;
        this.cacheMissBackup = Math.max(cacheSize / 1000, 1);
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     *
     * @param fasta The file to open.
     * @param index Pre-built FastaSequenceIndex, for the case in which one does not exist on disk.
     * @throws java.io.FileNotFoundException If the fasta or any of its supporting files cannot be found.
     */
    public CachingIndexedFastaSequenceFile(final File fasta, final FastaSequenceIndex index) {
        this(fasta, index, DEFAULT_CACHE_SIZE);
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     *
     * Looks for a index file for fasta on disk
     *
     * @param fasta The file to open.
     */
    public CachingIndexedFastaSequenceFile(final File fasta) throws FileNotFoundException {
        this(fasta, DEFAULT_CACHE_SIZE);
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     *
     * Looks for a index file for fasta on disk
     * Uses provided cacheSize instead of the default
     *
     * @param fasta The file to open.
     */
    public CachingIndexedFastaSequenceFile(final File fasta, final long cacheSize ) throws FileNotFoundException {
        super(fasta);
        if ( cacheSize < 0 ) throw new IllegalArgumentException("cacheSize must be > 0");
        this.cacheSize = cacheSize;
        this.cacheMissBackup = Math.max(cacheSize / 1000, 1);
    }

    /**
     * Print the efficiency (hits / queries) to logger with priority
     */
    public void printEfficiency(final Priority priority) {
        logger.log(priority, String.format("### CachingIndexedFastaReader: hits=%d misses=%d efficiency %.6f%%", cacheHits, cacheMisses, calcEfficiency()));
    }

    /**
     * Returns the efficiency (% of hits of all queries) of this object
     * @return
     */
    public double calcEfficiency() {
        return 100.0 * cacheHits / (cacheMisses + cacheHits * 1.0);
    }

    /**
     * @return the number of cache hits that have occurred
     */
    public long getCacheHits() {
        return cacheHits;
    }

    /**
     * @return the number of cache misses that have occurred
     */
    public long getCacheMisses() {
        return cacheMisses;
    }

    /**
     * @return the size of the cache we are using
     */
    public long getCacheSize() {
        return cacheSize;
    }

    /**
     * Gets the subsequence of the contig in the range [start,stop]
     *
     * Uses the sequence cache if possible, or updates the cache to handle the request.  If the range
     * is larger than the cache itself, just loads the sequence directly, not changing the cache at all
     *
     * @param contig Contig whose subsequence to retrieve.
     * @param start inclusive, 1-based start of region.
     * @param stop inclusive, 1-based stop of region.
     * @return The partial reference sequence associated with this range.
     */
    public ReferenceSequence getSubsequenceAt( final String contig, final long start, final long stop ) {
        final ReferenceSequence result;
        final Cache myCache = cache.get();

        if ( (stop - start) >= cacheSize ) {
            cacheMisses++;
            result = super.getSubsequenceAt(contig, start, stop);
        } else {
            // todo -- potential optimization is to check if contig.name == contig, as this in generally will be true
            SAMSequenceRecord contigInfo = super.getSequenceDictionary().getSequence(contig);

            if (stop > contigInfo.getSequenceLength())
                throw new PicardException("Query asks for data past end of contig");

            if ( start < myCache.start || stop > myCache.stop || myCache.seq == null || myCache.seq.getContigIndex() != contigInfo.getSequenceIndex() ) {
                cacheMisses++;
                myCache.start = Math.max(start - cacheMissBackup, 0);
                myCache.stop  = Math.min(start + cacheSize + cacheMissBackup, contigInfo.getSequenceLength());
                myCache.seq   = super.getSubsequenceAt(contig, myCache.start, myCache.stop);
                //System.out.printf("New cache at %s %d-%d%n", contig, cacheStart, cacheStop);
            } else {
                cacheHits++;
            }

            // at this point we determine where in the cache we want to extract the requested subsequence
            final int cacheOffsetStart = (int)(start - myCache.start);
            final int cacheOffsetStop = (int)(stop - start + cacheOffsetStart + 1);

            try {
                result = new ReferenceSequence(myCache.seq.getName(), myCache.seq.getContigIndex(), Arrays.copyOfRange(myCache.seq.getBases(), cacheOffsetStart, cacheOffsetStop));
            } catch ( ArrayIndexOutOfBoundsException e ) {
                throw new ReviewedStingException(String.format("BUG: bad array indexing.  Cache start %d and end %d, request start %d end %d, offset start %d and end %d, base size %d",
                        myCache.start, myCache.stop, start, stop, cacheOffsetStart, cacheOffsetStop, myCache.seq.getBases().length), e);
            }
        }

        if ( PRINT_EFFICIENCY && (getCacheHits() + getCacheMisses()) % PRINT_FREQUENCY == 0 )
            printEfficiency(Priority.INFO);
        return result;
    }
}