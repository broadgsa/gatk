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

package org.broadinstitute.gatk.engine.alignment.reference.bwt;

/**
 * Models a block of bases within the BWT.
 */
public class SequenceBlock {
    /**
     * Start position of this sequence within the BWT.
     */
    public final int sequenceStart;

    /**
     * Length of this sequence within the BWT.
     */
    public final int sequenceLength;


    /**
     * Occurrences of each letter up to this sequence block.
     */
    public final Counts occurrences;

    /**
     * Sequence for this segment.
     */
    public final byte[] sequence;

    /**
     * Create a new block within this BWT.
     * @param sequenceStart Starting position of this sequence within the BWT.
     * @param sequenceLength Length of this sequence.
     * @param occurrences How many of each base has been seen before this sequence began.
     * @param sequence The actual sequence from the BWT.
     */
    public SequenceBlock( int sequenceStart, int sequenceLength, Counts occurrences, byte[] sequence ) {
        this.sequenceStart = sequenceStart;
        this.sequenceLength = sequenceLength;
        this.occurrences = occurrences;
        this.sequence = sequence;
    }
}