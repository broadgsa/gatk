package net.sf.samtools;

import java.util.List;

/**
 * An individual bin in a BAM file, divided into chunks plus a linear index. 
 *
 * @author mhanna
 * @version 0.1
 */
public class Bin {
    /**
     * The reference sequence associated with this bin.
     */
    public final int referenceSequence;

    /**
     * The number of this bin within the BAM file.
     */
    public final int binNumber;

    /**
     * The chunks contained within this bin.
     */
    public final List<Chunk> chunks;

    public Bin(int referenceSequence, int binNumber, List<Chunk> chunks) {
        this.referenceSequence = referenceSequence;
        this.binNumber = binNumber;
        this.chunks = chunks;
    }
}
