package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.alignment.bwa.bwt.BWT;
import org.broadinstitute.sting.alignment.bwa.bwt.BWTReader;
import org.broadinstitute.sting.alignment.bwa.bwt.Base;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.alignment.Aligner;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.File;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Collections;
import java.util.EnumSet;

import net.sf.samtools.SAMRecord;
import net.sf.picard.util.SequenceUtil;

/**
 * Create imperfect alignments from the read to the genome represented by the given BWT / suffix array. 
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAAligner implements Aligner {
    /**
     * BWT in the forward direction.
     */
    private BWT forwardBWT;

    /**
     * BWT in the reverse direction.
     */
    private BWT reverseBWT;

    /**
     * Maximum edit distance (-n option from original BWA).
     */
    private static final int MAXIMUM_EDIT_DISTANCE = 4;

    /**
     * Maximum number of gap opens (-o option from original BWA).
     */
    private static final int MAXIMUM_GAP_OPENS = 1;

    /**
     * Maximum number of gap extensions (-e option from original BWA).
     */
    private static final int MAXIMUM_GAP_EXTENSIONS = -1;

    /**
     * Penalty for straight mismatches (-M option from original BWA).
     */
    private static final int MISMATCH_PENALTY = 3;

    /**
     * Penalty for gap opens (-O option from original BWA).
     */
    private static final int GAP_OPEN_PENALTY = 11;

    /**
     * Penalty for gap extensions (-E option from original BWA).
     */
    private static final int GAP_EXTENSION_PENALTY = 4;

    public BWAAligner( File forwardBWTFile, File reverseBWTFile ) {
        forwardBWT = new BWTReader(forwardBWTFile).read();
        reverseBWT = new BWTReader(reverseBWTFile).read();
    }

    public List<Alignment> align( SAMRecord read ) {
        byte[] bases = read.getReadBases();

        List<LowerBound> reverseLowerBounds = LowerBound.create(bases,reverseBWT);
//        for( int i = 0; i < reverseLowerBounds.size(); i++ )
//            System.out.printf("ReverseBWT: lb[%d] = %s%n",i,reverseLowerBounds.get(i));

        PriorityQueue<BWAAlignment> alignments = new PriorityQueue<BWAAlignment>();

        // Create a fictional initial alignment, with the position just off the end of the read, and the limits
        // set as the entire BWT.
        BWAAlignment initial = new BWAAlignment();
        initial.position = read.getReadLength();
        initial.loBound = 0;
        initial.hiBound = forwardBWT.length();
        initial.mismatches = 0;

        alignments.add(initial);

        while(!alignments.isEmpty()) {
            BWAAlignment alignment = alignments.remove();

            // Done with this particular alignment.
            if(alignment.position == 0)
                return Collections.<Alignment>singletonList(alignment);

            //System.out.printf("Processing alignments; queue size = %d, alignment = %s, bound = %d%n", alignments.size(), alignment, lowerBounds.get(alignment.position-1).value);            

            // if z < D(i) then return {}
            if( alignment.mismatches > reverseLowerBounds.get(alignment.position-1).value )
                continue;

            if( alignment.mismatches > MAXIMUM_EDIT_DISTANCE )
                continue;

            // For each base in { A, C, G, T }
            for(Base base: EnumSet.allOf(Base.class)) {
                // Create and initialize a new alignment, given that base as the candidate.
                BWAAlignment newAlignment = new BWAAlignment();
                newAlignment.position = alignment.position - 1;
                newAlignment.loBound = forwardBWT.counts(base) + forwardBWT.occurrences(base,alignment.loBound-1) + 1;
                newAlignment.hiBound = forwardBWT.counts(base) + forwardBWT.occurrences(base,alignment.hiBound);
                newAlignment.mismatches = alignment.mismatches;
                if( base.toASCII() != read.getReadBases()[newAlignment.position] )
                    newAlignment.mismatches++;

                // If this alignment is valid, add it to the list.
                if( newAlignment.loBound <= newAlignment.hiBound )
                    alignments.add(newAlignment);
            }
        }

        return Collections.emptyList();
    }

}
