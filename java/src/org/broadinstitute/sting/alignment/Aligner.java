package org.broadinstitute.sting.alignment;

import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * Create perfect alignments from the read to the genome represented by the given BWT / suffix array. 
 *
 * @author mhanna
 * @version 0.1
 */
public interface Aligner {
    /**
     * Align the read to the reference.
     * @param read Read to align.
     * @return A list of the alignments.
     */
    public List<Alignment> align(SAMRecord read);
}

