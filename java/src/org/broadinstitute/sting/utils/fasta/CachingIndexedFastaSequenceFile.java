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
import net.sf.picard.reference.*;
import net.sf.samtools.SAMSequenceRecord;

import java.io.File;
import java.util.Arrays;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * A caching version of the IndexedFastaSequenceFile that avoids going to disk as often as the raw indexer.
 *
 * Thread-safe!  Uses a lock object to protect write and access to the cache.
 */
public class CachingIndexedFastaSequenceFile extends IndexedFastaSequenceFile {
    private static final boolean USE_CACHE = true;
    private static final boolean PRINT_EFFICIENCY = false;
    private static final int PRINT_FREQUENCY = 10000;
    private static final long DEFAULT_CACHE_SIZE = 1000000;

    private static long ourStaticCacheSize = DEFAULT_CACHE_SIZE;

    long cacheHits = 0;
    long cacheMisses = 0;

    private static class Cache {
        long cacheStart = -1;
        long cacheStop = -1;
        long cacheSize = DEFAULT_CACHE_SIZE;
        long cacheMissBackup = 100;
        ReferenceSequence cachedSequence = null;

        public Cache(long cacheSize) {
            this.cacheSize = cacheSize;
            this.cacheMissBackup = Math.max(cacheSize / 100, 1);
        }
    }

    private static ThreadLocal<Cache> cache = new ThreadLocal<Cache> () {
        @Override protected Cache initialValue() {
            return new Cache(ourStaticCacheSize);
        }
    };

    /**
     * Same as general constructor but allows one to override the default cacheSize
     * @param file
     * @param index
     * @param cacheSize
     */
    public CachingIndexedFastaSequenceFile(final File file, final FastaSequenceIndex index, long cacheSize) {
        super(file, index);
        ourStaticCacheSize = cacheSize; // BUG: class only supports one univeral cache size
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     * @param file The file to open.
     * @param index Pre-built FastaSequenceIndex, for the case in which one does not exist on disk.
     * @throws java.io.FileNotFoundException If the fasta or any of its supporting files cannot be found.
     */
    public CachingIndexedFastaSequenceFile(final File file, final FastaSequenceIndex index) {
        this(file, index, DEFAULT_CACHE_SIZE);
    }

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     * @param file The file to open.
     */
    public CachingIndexedFastaSequenceFile(final File file) {
        this(file, DEFAULT_CACHE_SIZE);
    }

    public CachingIndexedFastaSequenceFile(final File file, long cacheSize ) {
        super(file);
        ourStaticCacheSize = cacheSize; // BUG: class only supports one univeral cache size
    }

    public void printEfficiency() {
        // comment out to disable tracking
        if ( (cacheHits + cacheMisses) % PRINT_FREQUENCY == 0 ) {
            System.out.printf("### CachingIndexedFastaReader: hits=%d misses=%d efficiency %.6f%%%n", cacheHits, cacheMisses, calcEfficiency());
        }
    }

    public double calcEfficiency() {
        return 100.0 * cacheHits / (cacheMisses + cacheHits * 1.0);
    }

    public long getCacheHits() {
        return cacheHits;
    }

    public long getCacheMisses() {
        return cacheMisses;
    }


    /**
     * Gets the subsequence of the contig in the range [start,stop]
     * @param contig Contig whose subsequence to retrieve.
     * @param start inclusive, 1-based start of region.
     * @param stop inclusive, 1-based stop of region.
     * @return The partial reference sequence associated with this range.
     */
    public ReferenceSequence getSubsequenceAt( String contig, long start, long stop ) {
        ReferenceSequence result;
        Cache myCache = cache.get();
        //System.out.printf("getSubsequentAt cache=%s%n", myCache);

        if ( ! USE_CACHE || (stop - start) >= myCache.cacheSize ) {
            cacheMisses++;
            result = super.getSubsequenceAt(contig, start, stop);
        } else {
            // todo -- potential optimization is to check if contig.name == contig, as this in generally will be true
            SAMSequenceRecord contigInfo = super.getSequenceDictionary().getSequence(contig);

            if (stop > contigInfo.getSequenceLength())
                throw new PicardException("Query asks for data past end of contig");

            if ( start < myCache.cacheStart || stop > myCache.cacheStop || myCache.cachedSequence == null || myCache.cachedSequence.getContigIndex() != contigInfo.getSequenceIndex() ) {
                cacheMisses++;
                myCache.cacheStart = Math.max(start - myCache.cacheMissBackup, 0);
                myCache.cacheStop = Math.min(myCache.cacheStart + myCache.cacheSize, contigInfo.getSequenceLength());
                myCache.cachedSequence = super.getSubsequenceAt(contig, myCache.cacheStart, myCache.cacheStop);
                //System.out.printf("New cache at %s %d-%d%n", contig, cacheStart, cacheStop);
            } else {
                cacheHits++;
            }

            // at this point we determine where in the cache we want to extract the requested subsequence
            int cacheOffsetStart = (int)(start - myCache.cacheStart);
            int cacheOffsetStop = (int)(stop - start + cacheOffsetStart + 1);

            try {
                result = new ReferenceSequence(myCache.cachedSequence.getName(), myCache.cachedSequence.getContigIndex(), Arrays.copyOfRange(myCache.cachedSequence.getBases(), cacheOffsetStart, cacheOffsetStop));
            } catch ( ArrayIndexOutOfBoundsException e ) {
                throw new ReviewedStingException(String.format("BUG: bad array indexing.  Cache start %d and end %d, request start %d end %d, offset start %d and end %d, base size %d",
                        myCache.cacheStart, myCache.cacheStop, start, stop, cacheOffsetStart, cacheOffsetStop, myCache.cachedSequence.getBases().length), e);
            }
        }

//        // comment out to disable testing
//        ReferenceSequence verify = super.getSubsequenceAt(contig, start, stop);
//        if ( ! Arrays.equals(verify.getBases(), result.getBases()) )
//            throw new ReviewedStingException(String.format("BUG: cached reference sequence not the same as clean fetched version at %s %d %d", contig, start, stop));

        if ( PRINT_EFFICIENCY ) printEfficiency();
        return result;
    }
}