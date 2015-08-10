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

package org.broadinstitute.gatk.utils.codecs.samread;

import htsjdk.tribble.Feature;

/**
 * Represents a SAM record read from a SAM text format file. 
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMReadFeature implements Feature {
    /**
     * Name of this read.
     */
    private final String readName;

    /**
     * Flags associated with this read.
     */
    private final int flags;

    /**
     * Contig to which this read is aligned.
     */
    private final String contig;

    /**
     * Position on contig to which this read is aligned.
     */
    private final int alignmentStart;

    /**
     * Position on contig at which this alignment ends.
     */
    private final int alignmentEnd;

    /**
     * Mapping quality for the read.
     */
    private final int mapQ;

    /**
     * Cigar string matching read to reference.
     */
    private final String cigarString;

    /**
     * Contig to which this read's pair is aligned.
     */
    private final String mateContig;

    /**
     * Position in contig to which this read's pair is aligned.
     */
    private final int mateAlignmentStart;

    /**
     * Size between pairs.
     */
    private final int insertSize;

    /**
     * Bases in this read.
     */
    private final byte[] bases;

    /**
     * Qualities constituting this read.
     */
    private final byte[] qualities;

    // Tags are not currently supported.

    /**
     * create the read feature.  Default protection so that only other classes in this package can create it.
     */
    SAMReadFeature(final String readName,
                   final int flags,
                   final String contig,
                   final int alignmentStart,
                   final int alignmentEnd,
                   final int mapQ,
                   final String cigarString,
                   final String mateContig,
                   final int mateAlignmentStart,
                   final int insertSize,
                   final byte[] bases,
                   final byte[] qualities) {
        this.readName = readName;
        this.flags = flags;
        this.contig = contig;
        this.alignmentStart = alignmentStart;
        this.alignmentEnd = alignmentEnd;
        this.mapQ = mapQ;
        this.cigarString = cigarString;
        this.mateContig = mateContig;
        this.mateAlignmentStart = mateAlignmentStart;
        this.insertSize = insertSize;
        this.bases = bases;
        this.qualities = qualities;
    }

    public String getReadName() {
        return readName;
    }

    public int getFlags() {
        return flags;
    }

    public String getReferenceName() {
        return contig;
    }

    public int getAlignmentStart() {
        return alignmentStart;
    }

    public int getAlignmentEnd() {
        return alignmentEnd;
    }

    /**
     * An alias for getReferenceName, required by Feature interface.
     * @return Aligned contig name.
     */
    public String getChr() {
        return getContig();
    }

    /**
     * An alias for getReferenceName, required by Feature interface.
     * @return Aligned contig name.
     */
    public String getContig() {
        return getReferenceName();
    }

    /**
     * An alias for getAlignmentEnd(), required by Feature interface.
     * @return End of alignment, inclusive.
     */
    public int getStart() {
        return getAlignmentStart();
    }

    /**
     * An alias for getAlignmentStart(), required by Feature interface.
     * @return Aligned position.  1-based.
     */
    public int getEnd() {
        return getAlignmentEnd();
    }    

    public int getMappingQuality() {
        return mapQ;
    }

    public String getCigarString() {
        return cigarString;
    }

    public String getMateReferenceName() {
        return mateContig;
    }

    public int getMateAlignmentStart() {
        return mateAlignmentStart;
    }

    public int getInferredInsertSize() {
        return insertSize;
    }

    public byte[] getReadBases() {
        return bases;    
    }

    public byte[] getReadQualities() {
        return qualities;
    }
}
