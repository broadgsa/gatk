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

package org.broadinstitute.gatk.engine.alignment;

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

/**
 * Represents an alignment of a read to a site in the reference genome.
 *
 * @author mhanna
 * @version 0.1
 */
public class Alignment {
    protected int contigIndex;
    protected long alignmentStart;
    protected boolean negativeStrand;
    protected int mappingQuality;

    protected char[] cigarOperators;
    protected int[] cigarLengths;

    protected int editDistance;
    protected String mismatchingPositions;

    protected int numMismatches;
    protected int numGapOpens;
    protected int numGapExtensions;
    protected int bestCount;
    protected int secondBestCount;

    /**
     * Gets the index of the given contig.
     * @return the inde
     */
    public int getContigIndex() { return contigIndex; }

    /**
     * Gets the starting position for the given alignment.
     * @return Starting position.
     */
    public long getAlignmentStart() { return alignmentStart; }

    /**
     * Is the given alignment on the reverse strand?
     * @return True if the alignment is on the reverse strand.
     */
    public boolean isNegativeStrand() { return negativeStrand; }

    /**
     * Gets the score of this alignment.
     * @return The score.
     */
    public int getMappingQuality() { return mappingQuality; }

    /**
     * Gets the edit distance; will eventually end up in the NM SAM tag
     * if this alignment makes it that far.
     * @return The edit distance.
     */
    public int getEditDistance() { return editDistance; }

    /**
     * A string representation of which positions mismatch; contents of MD tag.
     * @return String representation of mismatching positions.
     */
    public String getMismatchingPositions() { return mismatchingPositions; }
    
    /**
     * Gets the number of mismatches in the read.
     * @return Number of mismatches.
     */
    public int getNumMismatches() { return numMismatches; }

    /**
     * Get the number of gap opens.
     * @return Number of gap opens.
     */
    public int getNumGapOpens() { return numGapOpens; }

    /**
     * Get the number of gap extensions.
     * @return Number of gap extensions.
     */
    public int getNumGapExtensions() { return numGapExtensions; }

    /**
     * Get the number of best alignments.
     * @return Number of top scoring alignments.
     */
    public int getBestCount() { return bestCount; }

    /**
     * Get the number of second best alignments.
     * @return Number of second best scoring alignments.
     */
    public int getSecondBestCount() { return secondBestCount; }

    /**
     * Gets the cigar for this alignment.
     * @return sam-jdk formatted alignment.
     */
    public Cigar getCigar() {
        Cigar cigar = new Cigar();
        for(int i = 0; i < cigarOperators.length; i++) {
            CigarOperator operator = CigarOperator.characterToEnum(cigarOperators[i]);
            cigar.add(new CigarElement(cigarLengths[i],operator));
        }
        return cigar;
    }

    /**
     * Temporarily implement getCigarString() for debugging; the TextCigarCodec is unfortunately
     * package-protected.
     * @return
     */
    public String getCigarString() {
        Cigar cigar = getCigar();
        if(cigar.isEmpty()) return "*";

        StringBuilder cigarString = new StringBuilder();
        for(CigarElement element: cigar.getCigarElements()) {
            cigarString.append(element.getLength());
            cigarString.append(element.getOperator());
        }
        return cigarString.toString();
    }

    /**
     * Stub for inheritance.
     */
    public Alignment() {}    

    /**
     * Create a new alignment object.
     * @param contigIndex The contig to which this read aligned.
     * @param alignmentStart The point within the contig to which this read aligned.
     * @param negativeStrand Forward or reverse alignment of the given read.
     * @param mappingQuality How good does BWA think this mapping is?
     * @param cigarOperators The ordered operators in the cigar string.
     * @param cigarLengths The lengths to which each operator applies.
     * @param editDistance The edit distance (cumulative) of the read.
     * @param mismatchingPositions String representation of which bases in the read mismatch.
     * @param numMismatches Number of total mismatches in the read.
     * @param numGapOpens Number of gap opens in the read.
     * @param numGapExtensions Number of gap extensions in the read.
     * @param bestCount Number of best alignments in the read.
     * @param secondBestCount Number of second best alignments in the read.
     */
    public Alignment(int contigIndex,
                     int alignmentStart,
                     boolean negativeStrand,
                     int mappingQuality,
                     char[] cigarOperators,
                     int[] cigarLengths,
                     int editDistance,
                     String mismatchingPositions,
                     int numMismatches,
                     int numGapOpens,
                     int numGapExtensions,
                     int bestCount,
                     int secondBestCount) {
        this.contigIndex = contigIndex;
        this.alignmentStart = alignmentStart;
        this.negativeStrand = negativeStrand;
        this.mappingQuality = mappingQuality;
        this.cigarOperators = cigarOperators;
        this.cigarLengths = cigarLengths;
        this.editDistance = editDistance;
        this.mismatchingPositions = mismatchingPositions;
        this.numMismatches = numMismatches;
        this.numGapOpens = numGapOpens;
        this.numGapExtensions = numGapExtensions;
        this.bestCount = bestCount;
        this.secondBestCount = secondBestCount;
    }

    /**
     * Creates a read directly from an alignment.
     * @param alignment The alignment to convert to a read.
     * @param unmappedRead Source of the unmapped read.  Should have bases, quality scores, and flags.
     * @param newSAMHeader The new SAM header to use in creating this read.  Can be null, but if so, the sequence
     *                     dictionary in the
     * @return A mapped alignment.
     */
    public static SAMRecord convertToRead(Alignment alignment, SAMRecord unmappedRead, SAMFileHeader newSAMHeader) {
        SAMRecord read;
        try {
            read = (SAMRecord)unmappedRead.clone();
        }
        catch(CloneNotSupportedException ex) {
            throw new ReviewedGATKException("Unable to create aligned read from template.");
        }

        if(newSAMHeader != null)
            read.setHeader(newSAMHeader);

        // If we're realigning a previously aligned record, strip out the placement of the alignment.
        read.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        read.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        read.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        read.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);        

        if(alignment != null) {
            read.setReadUnmappedFlag(false);
            read.setReferenceIndex(alignment.getContigIndex());
            read.setAlignmentStart((int)alignment.getAlignmentStart());
            read.setReadNegativeStrandFlag(alignment.isNegativeStrand());
            read.setMappingQuality(alignment.getMappingQuality());
            read.setCigar(alignment.getCigar());
            if(alignment.isNegativeStrand()) {
                read.setReadBases(BaseUtils.simpleReverseComplement(read.getReadBases()));
                read.setBaseQualities(Utils.reverse(read.getBaseQualities()));
            }
            read.setAttribute("NM",alignment.getEditDistance());
            read.setAttribute("MD",alignment.getMismatchingPositions());
        }

        return read;
    }    
}
