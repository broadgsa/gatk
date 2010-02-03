package net.sf.samtools;

/**
 * The linear index associated with a given reference in a BAM index.
 *
 * @author mhanna
 * @version 0.1
 */
public class LinearIndex {
    /**
     * The reference sequence number for this linear index.
     */
    public final int referenceSequence;

    /**
     * The linear index entries within this bin.
     */
    public final long[] indexEntries;

    public LinearIndex(final int referenceSequence, final long[] indexEntries) {
        this.referenceSequence = referenceSequence;
        this.indexEntries = indexEntries;
    }
}
