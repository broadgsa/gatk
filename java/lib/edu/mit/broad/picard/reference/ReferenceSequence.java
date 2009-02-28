package edu.mit.broad.picard.reference;

/**
 * Wrapper around a reference sequence that has been read from a reference file.
 *
 * @author Tim Fennell
 */
public class ReferenceSequence {
    private String name;
    private byte[] bases;
    private int contigIndex;
    private int length;

    /**
     * Package level constructor that creates a fully formed ReferenceSequence
     *
     * @param name the name of the sequence from the source file
     * @param index the zero based index of this contig in the source file
     * @param bases the bases themselves stored as one-byte characters
     */
    ReferenceSequence(String name, int index, byte[] bases) {
        this.name = name;
        this.contigIndex = index;
        this.bases = bases;
        this.length = bases.length;
    }

    /** Gets the set of names given to this sequence in the source file. */
    public String getName() { return name; }

    /**
     * Gets the array of bases that define this sequence. The bases can include any
     * letter and possibly include masking information in the form of lower case
     * letters.  This array is mutable (obviously!) and it NOT a clone of the array
     * held interally.  Do not modify it!!!
     */
    public byte[] getBases() { return bases; }

    /** Gets the 0-based index of this contig in the source file from which it came. */
    public int getContigIndex() { return contigIndex; }

    /** Gets the length of this reference sequence in bases. */
    public int length() { return length; }
    
    public String toString() {
        return "ReferenceSequence " + getName();
    }
}
