package org.broadinstitute.sting.bwa;

/**
 * An alignment object to be used incrementally as the BWA aligner
 * inspects the read.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAAlignment implements Alignment {
    /**
     * Working variable.  The lower bound of the alignment within the BWT.
     */
    protected int loBound;

    /**
     * Working variable.  The upper bound of the alignment within the BWT.
     */
    protected int hiBound;

    /**
     * Current score for this alignment.
     */
    protected int score;

    /**
     * Gets the BWA score of this alignment.
     * @return BWA-style scores.  0 is best.
     */
    public int getScore() {
        return score;
    }

    /**
     * Compare this alignment to another alignment.
     * @param other Other alignment to which to compare.
     * @return < 0 if this < other, == 0 if this == other, > 0 if this > other
     */
    public int compareTo( Alignment other ) {
        return Integer.valueOf(score).compareTo(other.getScore());
    }
}
