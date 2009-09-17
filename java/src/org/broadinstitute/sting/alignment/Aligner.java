package org.broadinstitute.sting.alignment;

import net.sf.samtools.SAMRecord;

import java.io.File;
import java.util.List;

import org.broadinstitute.sting.alignment.bwa.*;
import org.broadinstitute.sting.alignment.bwa.bwt.BWT;
import org.broadinstitute.sting.alignment.bwa.bwt.SuffixArray;
import org.broadinstitute.sting.alignment.bwa.bwt.BWTReader;
import org.broadinstitute.sting.alignment.bwa.bwt.SuffixArrayReader;

/**
 * Create perfect alignments from the read to the genome represented by the given BWT / suffix array. 
 *
 * @author mhanna
 * @version 0.1
 */
public class Aligner {
    /**
     * BWT in the forward direction.
     */
    private BWT forwardBWT;

    /**
     * Suffix array in the forward direction.
     */
    private SuffixArray forwardSuffixArray;

    /**
     * BWT in the reverse direction.
     */
    private BWT reverseBWT;

    /**
     * Suffix array in the reverse direction.
     */
    private SuffixArray reverseSuffixArray;

    public Aligner( File forwardBWTFile, File forwardSuffixArrayFile, File reverseBWTFile, File reverseSuffixArrayFile ) {
        forwardBWT = new BWTReader(forwardBWTFile).read();
        forwardSuffixArray = new SuffixArrayReader(forwardSuffixArrayFile).read();

        reverseBWT = new BWTReader(reverseBWTFile).read();
        reverseSuffixArray = new SuffixArrayReader(reverseSuffixArrayFile).read();
    }

    /**
     * Align the read to the given reference.
     * @param read Read to align.
     * @return A list of the alignments.
     */
    public List<Alignment> align( SAMRecord read ) {
        List<LowerBound> lowerBounds = LowerBound.create(read,reverseBWT);
        return null;    
    }
}

