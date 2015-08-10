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

package org.broadinstitute.gatk.engine.alignment.bwa.java;

import htsjdk.samtools.Cigar;
import org.broadinstitute.gatk.engine.alignment.Alignment;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

/**
 * An alignment object to be used incrementally as the BWA aligner
 * inspects the read.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAAlignment extends Alignment implements Cloneable {
    /**
     * Track the number of alignments that have been created.
     */
    private static long numCreated;

    /**
     * Which number alignment is this?
     */
    private long creationNumber;

    /**
     * The aligner performing the alignments.
     */
    protected BWAJavaAligner aligner;

    /**
     * The sequence of matches/mismatches/insertions/deletions.
     */
    private AlignmentMatchSequence alignmentMatchSequence = new AlignmentMatchSequence();

    /**
     * Working variable.  How many bases have been matched at this point.
     */
    protected int position;

    /**
     * Working variable.  How many mismatches have been encountered at this point.
     */
    private int mismatches;

    /**
     * Number of gap opens in alignment.
     */
    private int gapOpens;

    /**
     * Number of gap extensions in alignment.
     */
    private int gapExtensions;

    /**
     * Working variable.  The lower bound of the alignment within the BWT.
     */
    protected long loBound;

    /**
     * Working variable.  The upper bound of the alignment within the BWT.
     */
    protected long hiBound;

    protected void setAlignmentStart(long position) {
        this.alignmentStart = position;
    }

    protected void setNegativeStrand(boolean negativeStrand) {
        this.negativeStrand = negativeStrand;
    }

    /**
     * Cache the score.
     */
    private int score;

    public Cigar getCigar() {
        return alignmentMatchSequence.convertToCigar(isNegativeStrand());
    }

    /**
     * Gets the current state of this alignment (state of the last base viewed)..
     * @return Current state of the alignment.
     */
    public AlignmentState getCurrentState() {
        return alignmentMatchSequence.getCurrentState();
    }

    /**
     * Adds the given state to the current alignment.
     * @param state State to add to the given alignment.
     */
    public void addState( AlignmentState state ) {
        alignmentMatchSequence.addNext(state);    
    }

    /**
     * Gets the BWA score of this alignment.
     * @return BWA-style scores.  0 is best.
     */
    public int getScore() {
        return score;
    }

    public int getMismatches() { return mismatches; }
    public int getGapOpens() { return gapOpens; }
    public int getGapExtensions() { return gapExtensions; }

    public void incrementMismatches() {
        this.mismatches++;
        updateScore();
    }

    public void incrementGapOpens() {
        this.gapOpens++;
        updateScore();
    }

    public void incrementGapExtensions() {
        this.gapExtensions++;
        updateScore();
    }

    /**
     * Updates the score based on new information about matches / mismatches.
     */
    private void updateScore() {
        score = mismatches*aligner.MISMATCH_PENALTY + gapOpens*aligner.GAP_OPEN_PENALTY + gapExtensions*aligner.GAP_EXTENSION_PENALTY;
    }

    /**
     * Create a new alignment with the given parent aligner.
     * @param aligner Aligner being used.
     */
    public BWAAlignment( BWAJavaAligner aligner ) {
        this.aligner = aligner;
        this.creationNumber = numCreated++;
    }

    /**
     * Clone the alignment.
     * @return New instance of the alignment.
     */
    public BWAAlignment clone() {
        BWAAlignment newAlignment = null;
        try {
            newAlignment = (BWAAlignment)super.clone();
        }
        catch( CloneNotSupportedException ex ) {
            throw new ReviewedGATKException("Unable to clone BWAAlignment.");
        }
        newAlignment.creationNumber = numCreated++;
        newAlignment.alignmentMatchSequence = alignmentMatchSequence.clone();

        return newAlignment;
    }

    /**
     * How many bases in the read match the given state.
     * @param state State to test.
     * @return number of bases which match that state.
     */
    public int getNumberOfBasesMatchingState(AlignmentState state) {
        return alignmentMatchSequence.getNumberOfBasesMatchingState(state);
    }

    /**
     * Compare this alignment to another alignment.
     * @param rhs Other alignment to which to compare.
     * @return < 0 if this < other, == 0 if this == other, > 0 if this > other
     */
    public int compareTo(Alignment rhs) {
        BWAAlignment other = (BWAAlignment)rhs;

        // If the scores are different, disambiguate using the score.
        if(score != other.score)
            return score > other.score ? 1 : -1;

        // Otherwise, use the order in which the elements were created.
        if(creationNumber != other.creationNumber)
            return creationNumber > other.creationNumber ? -1 : 1;

        return 0;
    }

    public String toString() {
        return String.format("position: %d, strand: %b, state: %s, mismatches: %d, gap opens: %d, gap extensions: %d, loBound: %d, hiBound: %d, score: %d, creationNumber: %d", position, negativeStrand, alignmentMatchSequence.getCurrentState(), mismatches, gapOpens, gapExtensions, loBound, hiBound, getScore(), creationNumber);
    }
}
