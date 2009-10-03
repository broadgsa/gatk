package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.utils.StingException;

/**
 * An alignment object to be used incrementally as the BWA aligner
 * inspects the read.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAAlignment implements Alignment, Cloneable {
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
    protected BWAAligner aligner;

    /**
     * Start of the final alignment.
     */
    protected int alignmentStart;

    /**
     * Is this match being treated as a negative or positive strand?
     */
    protected boolean negativeStrand;

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
    protected int mismatches;

    /**
     * Number of gap opens in alignment.
     */
    protected int gapOpens;

    /**
     * Number of gap extensions in alignment.
     */
    protected int gapExtensions;

    /**
     * Working variable.  The lower bound of the alignment within the BWT.
     */
    protected int loBound;

    /**
     * Working variable.  The upper bound of the alignment within the BWT.
     */
    protected int hiBound;

    /**
     * Gets the starting position for the given alignment.
     * @return Starting position.
     */
    public int getAlignmentStart() {
        return alignmentStart;
    }

    /**
     * Is the given alignment on the reverse strand?
     * @return True if the alignment is on the reverse strand.
     */
    public boolean isNegativeStrand() {
        return negativeStrand;    
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
        return mismatches*aligner.MISMATCH_PENALTY + gapOpens*aligner.GAP_OPEN_PENALTY + gapExtensions*aligner.GAP_EXTENSION_PENALTY;
    }

    public int getMismatches() { return mismatches; }
    public int getGapOpens() { return gapOpens; }
    public int getGapExtensions() { return gapExtensions; }

    /**
     * Create a new alignment with the given parent aligner.
     * @param aligner Aligner being used.
     */
    public BWAAlignment( BWAAligner aligner ) {
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
            throw new StingException("Unable to clone BWAAlignment.");
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

        // If the scores are equal, use the score to disambiguate.
        int scoreComparison = Integer.valueOf(getScore()).compareTo(other.getScore());
        if( scoreComparison != 0 )
            return scoreComparison;

        return -Long.valueOf(this.creationNumber).compareTo(other.creationNumber);
    }

    public String toString() {
        return String.format("position: %d, state: %s, mismatches: %d, gap opens: %d, gap extensions: %d, loBound: %d, hiBound: %d, score: %d, creationNumber: %d", position, alignmentMatchSequence.getCurrentState(), mismatches, gapOpens, gapExtensions, loBound, hiBound, getScore(), creationNumber);
    }
}
