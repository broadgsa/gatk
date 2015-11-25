/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;


/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *
 *  Class ArtificialPatternedSAMIterator
 *
 * This class allows you to pattern the artificial sam iterator, asking for reads
 * in order or out of order.
 */
public class ArtificialPatternedSAMIterator extends ArtificialSAMIterator {

    /** the pattern we're implementing */
    public enum PATTERN {
        RANDOM_READS, IN_ORDER_READS;
    }

    // our pattern
    private final PATTERN mPattern;

    /**
     * this is pretty heavy (and it could be extremely heavy, given the amount of reads they request, but it
     * allows us to give them each read once, reguardless of the order specified
     */
    private final int[] reads;
    private final int readCount;

    /**
     * create the fake iterator, given the mapping of chromosomes and read counts.  If pattern
     * is specified to be random, it will generate reads that are randomly placed on the current chromosome
     *
     * @param startingChr the starting chromosome
     * @param endingChr   the ending chromosome
     * @param readCount   the number of reads in each chromosome
     * @param header      the associated header
     * @param pattern     the pattern to implement
     */
    public ArtificialPatternedSAMIterator( int startingChr, int endingChr, int readCount, int unmappedReadCount, SAMFileHeader header, PATTERN pattern ) {
        super(startingChr, endingChr, readCount, unmappedReadCount, header);
        mPattern = pattern;
        this.readCount = readCount;
        reads = new int[readCount];

        for (int x = 0; x < readCount; x++) {
            reads[x] = x+1;
        }
        if (pattern == PATTERN.RANDOM_READS) {
            // scramble a bunch of the reads
            for (int y = 0; y < readCount; y++) {
                int ranOne = (int) Math.round(Math.random() * ( readCount - 1 ));
                int ranTwo = (int) Math.round(Math.random() * ( readCount - 1 ));
                int temp = reads[ranOne];
                reads[ranOne] = reads[ranTwo];
                reads[ranTwo] = temp;
            }
            /**
             *  up to this point there's no garauntee that the random() has made the reads out of order (though it's
             *  extremely extremely unlikely it's failed).  Let's make sure there at least out of order:
             */
            if (this.reads[0] < this.reads[reads.length - 1]) {
                int temp = reads[0];
                reads[0] = reads[reads.length - 1];
                reads[reads.length - 1] = temp;
            }

        }

    }

    /**
     * override the default ArtificialSAMIterator createNextRead method, which creates the next read
     *
     * @return
     */
    protected boolean createNextRead() {
        if (currentRead > rCount) {
            currentChromo++;
            currentRead = 1;
        }
        // check for end condition, have we finished the chromosome listing, and have no unmapped reads
        if (currentChromo >= eChromosomeCount) {
            if (unmappedRemaining < 1) {
                this.next = null;
                return false;
            } else {
                ++totalReadCount;
                this.next = ArtificialSAMUtils.createArtificialRead(this.header,
                        String.valueOf(totalReadCount),
                        SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX,
                        SAMRecord.NO_ALIGNMENT_START,
                        50);
                --unmappedRemaining;
                return true;
            }
        }
        ++totalReadCount;
        this.next = getNextRecord(currentRead);

        ++currentRead;
        return true;
    }


    /**
     * get the next read, given it's index in the chromosome
     *
     * @param read the read index in the chromosome
     *
     * @return a SAMRecord
     */
    private SAMRecord getNextRecord( int read ) {
        if (read > this.readCount) {
            return ArtificialSAMUtils.createArtificialRead(this.header, String.valueOf(reads[readCount - 1]), currentChromo, reads[readCount - 1], 50);
        }
        return ArtificialSAMUtils.createArtificialRead(this.header, String.valueOf(reads[read-1]), currentChromo, reads[read-1], 50);
    }

}
