package org.broadinstitute.sting.alignment;

/**
 * Represents an alignment of a read to a site in the reference genome.
 *
 * @author mhanna
 * @version 0.1
 */
public interface Alignment extends Comparable<Alignment> {
    /**
     * Is the given alignment on the reverse strand?
     * @return True if the alignment is on the reverse strand.
     */
    public boolean isNegativeStrand();

    /**
     * Gets the starting position for the given alignment.
     * @return Starting position.
     */
    public long getAlignmentStart();

    /**
     * Gets the score of this alignment.
     * @return The score.
     */
    public int getScore();
}
