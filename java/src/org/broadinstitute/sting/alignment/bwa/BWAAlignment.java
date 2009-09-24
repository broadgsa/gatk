package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.alignment.Alignment;

/**
 * An alignment object to be used incrementally as the BWA aligner
 * inspects the read.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAAlignment implements Alignment {
    /**
     * Start of the final alignment.
     */
    protected int alignmentStart;

    /**
     * Working variable.  Is this match being treated as a negative or positive strand?
     */
    protected boolean negativeStrand;

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
     * Indicates the current state of an alignment.  Are we in an insertion?  Deletion?
     */
    protected AlignmentState state;

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
     * Gets the BWA score of this alignment.
     * @return BWA-style scores.  0 is best.
     */
    public int getScore() {
        return mismatches + gapOpens + gapExtensions;
    }

    public int getMismatches() { return mismatches; }
    public int getGapOpens() { return gapOpens; }
    public int getGapExtensions() { return gapExtensions; }   

    /**
     * Compare this alignment to another alignment.
     * @param other Other alignment to which to compare.
     * @return < 0 if this < other, == 0 if this == other, > 0 if this > other
     */
    public int compareTo(Alignment other) {
        // If the scores are equal, use the position to disambiguate order.
        int scoreComparison = Integer.valueOf(getScore()).compareTo(other.getScore());
        if( scoreComparison != 0 )
            return scoreComparison;
        else
            return -Integer.valueOf(position).compareTo(((BWAAlignment)other).position);
    }

    public String toString() {
        return String.format("position: %d, state: %s, mismatches: %d, gap opens: %d, gap extensions: %d, loBound: %d, hiBound: %d, score: %d", position, state, mismatches, gapOpens, gapExtensions, loBound, hiBound, getScore());
    }
}
