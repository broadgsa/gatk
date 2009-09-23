package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.alignment.bwa.bwt.*;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.alignment.Aligner;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.File;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Collections;
import java.util.EnumSet;

import net.sf.samtools.SAMRecord;

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
     * Suffix array in the forward direction.
     */
    private SuffixArray reverseSuffixArray;

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

    public BWAAligner( File forwardBWTFile, File reverseBWTFile, File reverseSuffixArrayFile ) {
        forwardBWT = new BWTReader(forwardBWTFile).read();
        reverseBWT = new BWTReader(reverseBWTFile).read();
        reverseSuffixArray = new SuffixArrayReader(reverseSuffixArrayFile).read();
    }

    public List<Alignment> align( SAMRecord read ) {
        byte[] forwardBases = read.getReadBases();
        byte[] reverseBases = BaseUtils.simpleReverseComplement(forwardBases);

        List<LowerBound> forwardLowerBounds = LowerBound.create(forwardBases,forwardBWT);
        //for( int i = 0; i < forwardLowerBounds.size(); i++ )
        //    System.out.printf("ForwardBWT: lb[%d] = %s%n",i,forwardLowerBounds.get(i));
        List<LowerBound> reverseLowerBounds = LowerBound.create(reverseBases,reverseBWT);
        //for( int i = 0; i < reverseLowerBounds.size(); i++ )
        //    System.out.printf("ReverseBWT: lb[%d] = %s%n",i,reverseLowerBounds.get(i));

        PriorityQueue<BWAAlignment> alignments = new PriorityQueue<BWAAlignment>();

        // Create a fictional initial alignment, with the position just off the end of the read, and the limits
        // set as the entire BWT.
        BWAAlignment initial = new BWAAlignment();
        initial.negativeStrand = true;
        initial.position = 0;
        initial.loBound = 0;
        initial.hiBound = forwardBWT.length();
        initial.mismatches = 0;

        BWAAlignment initialReverse = new BWAAlignment();
        initialReverse.negativeStrand = false;
        initialReverse.position = 0;
        initialReverse.loBound = 0;
        initialReverse.hiBound = reverseBWT.length();
        initialReverse.mismatches = 0;

        alignments.add(initial);
        //alignments.add(initialReverse);

        while(!alignments.isEmpty()) {
            BWAAlignment alignment = alignments.remove();

            byte[] bases = alignment.negativeStrand ? forwardBases : reverseBases;
            BWT bwt = alignment.negativeStrand ? reverseBWT : forwardBWT;
            List<LowerBound> lowerBounds = alignment.negativeStrand ? forwardLowerBounds : reverseLowerBounds;

            // Done with this particular alignment.
            if(alignment.position == read.getReadLength()-1) {
                alignment.alignmentStart = reverseBWT.length() - (reverseSuffixArray.get(alignment.loBound)+read.getReadLength()) + 1;
                return Collections.<Alignment>singletonList(alignment);
            }

            //System.out.printf("Processing alignments; queue size = %d, alignment = %s, bound = %d%n", alignments.size(), alignment, lowerBounds.get(alignment.position).value);

            // if z < D(i) then return {}
            int mismatches = MAXIMUM_EDIT_DISTANCE - alignment.mismatches;            
            if( mismatches < lowerBounds.get(alignment.position+1).value )
                continue;

            if( alignment.mismatches > MAXIMUM_EDIT_DISTANCE )
                continue;

            // For each base in { A, C, G, T }
            for(Base base: EnumSet.allOf(Base.class)) {
                // Create and initialize a new alignment, given that base as the candidate.
                BWAAlignment newAlignment = new BWAAlignment();
                newAlignment.negativeStrand = alignment.negativeStrand;
                newAlignment.position = alignment.position + 1;

                newAlignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
                newAlignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);
                newAlignment.mismatches = alignment.mismatches;
                if( base.toASCII() != bases[newAlignment.position] )
                    newAlignment.mismatches++;

                // If this alignment is valid, add it to the list.
                if( newAlignment.loBound <= newAlignment.hiBound )
                    alignments.add(newAlignment);
            }
        }

        return Collections.emptyList();
    }

}
