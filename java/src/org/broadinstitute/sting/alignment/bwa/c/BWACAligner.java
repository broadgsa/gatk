package org.broadinstitute.sting.alignment.bwa.c;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.alignment.Alignment;

import java.util.*;

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
     * Destroy the BWA/C thunk.
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
     * Get a iterator of alignments, batched by mapping quality.
     * @param bases List of bases.
     * @return Iterator to alignments.
     */
    public Iterator<Alignment[]> getAlignments(final byte[] bases) {
        final BWAPath[] paths = getPaths(bases);
        return new Iterator<Alignment[]>() {
            /**
             * The last position accessed.
             */
            private int position = 0;

            /**
             * Whether all alignments have been seen based on the current position.
             * @return True if any more alignments are pending.  False otherwise.
             */
            public boolean hasNext() { return position < paths.length; }

            /**
             * Return the next cross-section of alignments, based on mapping quality.
             * @return Array of the next set of alignments of a given mapping quality.
             */
            public Alignment[] next() {
                if(position >= paths.length)
                    throw new UnsupportedOperationException("Out of alignments to return.");
                int score = paths[position].score;
                int startingPosition = position;
                while(position < paths.length && paths[position].score == score) position++;
                return convertPathsToAlignments(bases,Arrays.copyOfRange(paths,startingPosition,position));
            }

            /**
             * Unsupported.
             */
            public void remove() { throw new UnsupportedOperationException("Cannot remove from an alignment iterator"); }
        };
    }

    /**
     * Get a iterator of aligned reads, batched by mapping quality.
     * @param read Read to align.
     * @return Iterator to alignments.
     */
    public Iterator<SAMRecord[]> align(final SAMRecord read) {
        final Iterator<Alignment[]> alignments = getAlignments(read.getReadBases());
        return new Iterator<SAMRecord[]>() {
            /**
             * Whether all alignments have been seen based on the current position.
             * @return True if any more alignments are pending.  False otherwise.
             */
            public boolean hasNext() { return alignments.hasNext(); }

            /**
             * Return the next cross-section of alignments, based on mapping quality.
             * @return Array of the next set of alignments of a given mapping quality.
             */
            public SAMRecord[] next() {
                Alignment[] alignmentsOfQuality = alignments.next();
                SAMRecord[] reads = new SAMRecord[alignmentsOfQuality.length];
                for(int i = 0; i < alignmentsOfQuality.length; i++) {
                    reads[i] = convertAlignmentToRead(read,alignmentsOfQuality[i]);    
                }
                return reads;
            }

            /**
             * Unsupported.
             */
            public void remove() { throw new UnsupportedOperationException("Cannot remove from an alignment iterator"); }
        };
    }

    /**
     * Push all alignment data into individual SAMRecords, gaining in convenience but losing some of
     * the additional data stored in an alignment object.
     * @param read The read to align.
     * @return A list of potential alignments.
     */
    public SAMRecord[] alignAll(SAMRecord read) {
        Alignment[] alignments = getAllAlignments(read.getReadBases());
        SAMRecord[] reads = new SAMRecord[alignments.length];
        for(int i = 0; i < alignments.length; i++)
            reads[i] = convertAlignmentToRead(read,alignments[i]);

        return reads;
    }

    /**
     * Creates a read directly from an alignment.
     * @param unmappedRead Source of the unmapped read.  Should have bases, quality scores, and flags.
     * @param alignment The target alignment for this read.
     * @return A mapped alignment.
     */
    public SAMRecord convertAlignmentToRead(SAMRecord unmappedRead, Alignment alignment) {
        SAMRecord read = null;
        try {
            read = (SAMRecord)unmappedRead.clone();
            read.setReadUmappedFlag(false);
            read.setAlignmentStart((int)alignment.getAlignmentStart());
            read.setReadNegativeStrandFlag(alignment.isNegativeStrand());
            read.setMappingQuality(alignment.getMappingQuality());
            read.setCigar(alignment.getCigar());
            if(alignment.isNegativeStrand()) {
                read.setReadBases(BaseUtils.reverse(read.getReadBases()));
                read.setBaseQualities(BaseUtils.reverse(read.getBaseQualities()));
            }
        }
        catch(CloneNotSupportedException ex) {
            throw new StingException("Unable to create aligned read from template.");
        }
        return read;
    }

    /**
     * Align the given base array to the BWT.  The base array should be in ASCII format.
     * @param bases ASCII representation of byte array.
     * @return an array of indices into the bwa.
     */
    public Alignment[] getAllAlignments(byte[] bases) {
        if(thunkPointer == 0)
            throw new StingException("BWA/C align attempted, but BWA/C is not properly initialized.");
        BWAPath[] paths = getPaths(bases);
        return convertPathsToAlignments(thunkPointer,bases,paths);
    }

    /**
     * Get the paths associated with the given base string.
     * @param bases List of bases.
     * @return A set of paths through the BWA.
     */
    public BWAPath[] getPaths(byte[] bases) {
        if(thunkPointer == 0)
            throw new StingException("BWA/C getPaths attempted, but BWA/C is not properly initialized.");
        BWAPath[] paths = getPaths(thunkPointer,bases);
        return paths;
    }

    /**
     * Do the extra steps involved in converting a local alignment to a global alignment.
     * @param bases ASCII representation of byte array.
     * @param paths Paths through the current BWT.
     * @return A list of alignments.
     */
    protected Alignment[] convertPathsToAlignments(byte[] bases, BWAPath[] paths) {
        if(thunkPointer == 0)
            throw new StingException("BWA/C convertPathsToAlignments attempted, but BWA/C is not properly initialized.");
        return convertPathsToAlignments(thunkPointer,bases,paths);
    }

    /**
     * Caller to the path generation functionality within BWA/C.  Call this method's getPaths() wrapper (above) instead.
     * @param thunkPointer pointer to the C++ object managing BWA/C.
     * @param bases ASCII representation of byte array.
     * @return A list of paths through the specified BWT.
     */
    protected native BWAPath[] getPaths(long thunkPointer, byte[] bases);

    /**
     * Do the extra steps involved in converting a local alignment to a global alignment.
     * Call this method's convertPathsToAlignments() wrapper (above) instead.
     * @param thunkPointer pointer to the C++ object managing BWA/C.
     * @param bases ASCII representation of byte array.
     * @param paths Paths through the current BWT.
     * @return A list of alignments.
     */
    protected native Alignment[] convertPathsToAlignments(long thunkPointer, byte[] bases, BWAPath[] paths);
}
