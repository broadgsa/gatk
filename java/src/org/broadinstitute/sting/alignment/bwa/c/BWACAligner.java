package org.broadinstitute.sting.alignment.bwa.c;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.alignment.bwa.c.BWACConfiguration;

/**
 * An aligner using the BWA/C implementation.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWACAligner {
    static {
        System.loadLibrary("bwa");
    }

    /**
     * A pointer to the C++ object representing the BWA engine.
     */
    private long thunkPointer = 0;

    /**
     * Create a pointer to the BWA/C thunk.
     * @param configuration Configuration of the aligner.
     * @return Pointer to the BWA/C thunk.
     */
    protected native long create(BWACConfiguration configuration);

    /**
     * Destroy the 
     * @param thunkPointer Pointer to the allocated thunk.
     */
    protected native void destroy(long thunkPointer);

    public BWACAligner(BWACConfiguration configuration) {
        if(thunkPointer != 0)
            throw new StingException("BWA/C attempting to reinitialize.");
        thunkPointer = create(configuration);
    }

    /**
     * Close this instance of the BWA pointer and delete its resources.
     */
    public void close() {
        if(thunkPointer == 0)
            throw new StingException("BWA/C close attempted, but BWA/C was never properly initialized.");
        destroy(thunkPointer);
    }

    /**
     * Align the given base array to the BWT.  The base array should be in ASCII format.
     * @param bases ASCII representation of byte array.
     * @return an array of indices into the bwa.
     */
    public Alignment[] getAlignments(byte[] bases) {
        if(thunkPointer == 0)
            throw new StingException("BWA/C align attempted, but BWA/C was never properly initialized.");
        return getAlignments(thunkPointer,bases);
    }

    /**
     * Push all alignment data into individual SAMRecords, gaining in convenience but losing some of
     * the additional data stored in an alignment object.
     * @param read The read to align.
     * @return A list of potential alignments.
     */
    public SAMRecord[] align(SAMRecord read) {
        Alignment[] alignments = getAlignments(read.getReadBases());
        SAMRecord[] reads = new SAMRecord[alignments.length];
        for(int i = 0; i < alignments.length; i++) {
            try {
                reads[i] = (SAMRecord)read.clone();
                reads[i].setReadUmappedFlag(false);
                reads[i].setAlignmentStart((int)alignments[i].getAlignmentStart());
                reads[i].setReadNegativeStrandFlag(alignments[i].isNegativeStrand());
                reads[i].setMappingQuality(alignments[i].getMappingQuality());
                reads[i].setCigar(alignments[i].getCigar());
                if(alignments[i].isNegativeStrand()) {
                    reads[i].setReadBases(BaseUtils.reverse(reads[i].getReadBases()));
                    reads[i].setBaseQualities(BaseUtils.reverse(reads[i].getBaseQualities()));
                }
            }
            catch(CloneNotSupportedException ex) {
                throw new StingException("Unable to create aligned read from template.");
            }
        }

        return reads;
    }

    /**
     * Caller to the .  The base array should be in
     * ASCII format.
     * @param thunkPointer pointer to the C++ object managing BWA/C.
     * @param bases ASCII representation of byte array.
     */
    protected native Alignment[] getAlignments(long thunkPointer, byte[] bases);
}
