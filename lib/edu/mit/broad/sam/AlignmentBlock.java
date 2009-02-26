package edu.mit.broad.sam;

/**
 * Represents the contiguous alignment of a subset of read bases to a reference
 * sequence. Simply put an alignment block tells you that read bases from
 * readStart are aligned to the reference (matching or mismatching) from
 * referenceStart for length bases.
 *
 * @author Tim Fennell
 */
public class AlignmentBlock {
    private int readStart;
    private int referenceStart;
    private int length;

    /** Constructs a new alignment block with the supplie read and ref starts and length. */
    AlignmentBlock(int readStart, int referenceStart, int length) {
        this.readStart = readStart;
        this.referenceStart = referenceStart;
        this.length = length;
    }

    /** The first, 1-based, base in the read that is aligned to the reference reference. */
    public int getReadStart() { return readStart; }

    /** The first, 1-based, position in the reference to which the read is aligned. */
    public int getReferenceStart() { return referenceStart; }

    /** The number of contiguous bases aligned to the reference. */
    public int getLength() { return length; }
}
