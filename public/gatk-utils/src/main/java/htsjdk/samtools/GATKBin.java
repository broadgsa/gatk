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

package htsjdk.samtools;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A temporary solution to work around Java access rights issues:
 * override GATKBin and make it public.
 * TODO: Eliminate once we determine the final fate of the BAM index reading code.
 */
public class GATKBin implements Comparable<GATKBin> {
    /**
     * The reference sequence associated with this bin.
     */
    private final int referenceSequence;

    /**
     * The number of this bin within the BAM file.
     */
    private final int binNumber;

    /**
     * The chunks associated with this bin.
     */
    private GATKChunk[] chunkList;

    public GATKBin(Bin bin) {
        this(bin.getReferenceSequence(),bin.getBinNumber());
    }

    public GATKBin(final int referenceSequence, final int binNumber) {
        this.referenceSequence = referenceSequence;
        this.binNumber = binNumber;
    }

    public int getReferenceSequence() {
        return referenceSequence;
    }

    public int getBinNumber() {
        return binNumber;
    }

    /**
     * Convert this GATKBin to a normal bin, for processing with the standard BAM query interface.
     * @return
     */
    public Bin toBin() {
        return new Bin(referenceSequence,binNumber);
    }

    /**
     * See whether two bins are equal.  If the ref seq and the bin number
     * are equal, assume equality of the chunk list.
     * @param other The other Bin to which to compare this.
     * @return True if the two bins are equal.  False otherwise.
     */
    @Override
    public boolean equals(Object other) {
        if(other == null) return false;
        if(!(other instanceof GATKBin)) return false;

        GATKBin otherBin = (GATKBin)other;
        return this.referenceSequence == otherBin.referenceSequence && this.binNumber == otherBin.binNumber;
    }

    /**
     * Compute a unique hash code for the given reference sequence and bin number.
     * @return A unique hash code.
     */
    @Override
    public int hashCode() {
        return ((Integer)referenceSequence).hashCode() ^ ((Integer)binNumber).hashCode();
    }

    /**
     * Compare two bins to see what ordering they should appear in.
     * @param other Other bin to which this bin should be compared.
     * @return -1 if this < other, 0 if this == other, 1 if this > other.
     */
    public int compareTo(GATKBin other) {
        if(other == null)
            throw new ClassCastException("Cannot compare to a null object");

        // Check the reference sequences first.
        if(this.referenceSequence != other.referenceSequence)
            return referenceSequence - other.referenceSequence;

        // Then check the bin ordering.
        return binNumber - other.binNumber;
    }

    /**
     * Sets the chunks associated with this bin
     */
    public void setChunkList(GATKChunk[] list){
        chunkList = list;
    }

    /**
     * Gets the list of chunks associated with this bin.
     * @return the chunks in this bin.  If no chunks are associated, an empty list will be returned.
     */
    public GATKChunk[] getChunkList(){
        if(chunkList == null)
            return new GATKChunk[0];
        return chunkList;
    }

    // HACK: Using this classes package permissions to further hack the CRAM created SAMRecord's indexing bin and binary attributes.
    public static Integer getReadIndexingBin(final SAMRecord read) {
        return read.getIndexingBin();
    }
    public static void setReadIndexingBin(final SAMRecord read, final Integer indexingBin) {
        read.setIndexingBin(indexingBin);
    }
    public static SAMBinaryTagAndValue getReadBinaryAttributes(final SAMRecord read) {
        return read.getBinaryAttributes();
    }
}
