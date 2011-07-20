package org.broadinstitute.sting.alignment;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * Create perfect alignments from the read to the genome represented by the given BWT / suffix array. 
 *
 * @author mhanna
 * @version 0.1
 */
public interface Aligner {
    /**
     * Close this instance of the BWA pointer and delete its resources.
     */
    public void close();    

    /**
     * Allow the aligner to choose one alignment randomly from the pile of best alignments.
     * @param bases Bases to align.
     * @return An align
     */
    public Alignment getBestAlignment(final byte[] bases);

    /**
     * Align the read to the reference.
     * @param read Read to align.
     * @param header Optional header to drop in place.
     * @return A list of the alignments.
     */
    public SAMRecord align(final SAMRecord read, final SAMFileHeader header);

    /**
     * Get a iterator of alignments, batched by mapping quality.
     * @param bases List of bases.
     * @return Iterator to alignments.
     */
    public Iterable<Alignment[]> getAllAlignments(final byte[] bases);

    /**
     * Get a iterator of aligned reads, batched by mapping quality.
     * @param read Read to align.
     * @param newHeader Optional new header to use when aligning the read.  If present, it must be null.
     * @return Iterator to alignments.
     */
    public Iterable<SAMRecord[]> alignAll(final SAMRecord read, final SAMFileHeader newHeader);
}


