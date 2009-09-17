package org.broadinstitute.sting.alignment;

/**
 * Represents an alignment of a read to a site in the reference genome.
 *
 * @author mhanna
 * @version 0.1
 */
public interface Alignment extends Comparable<Alignment> {
    /**
     * Gets the score of this alignment.
     * @return The score.
     */
    public int getScore();
}
