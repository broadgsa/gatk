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
    private static final int MAXIMUM_GAP_EXTENSIONS = 6;

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
        initial.state = AlignmentState.MATCH_MISMATCH;
        initial.mismatches = 0;

        BWAAlignment initialReverse = new BWAAlignment();
        initialReverse.negativeStrand = false;
        initialReverse.position = 0;
        initialReverse.loBound = 0;
        initialReverse.hiBound = reverseBWT.length();
        initialReverse.state = AlignmentState.MATCH_MISMATCH;
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
                alignment.alignmentStart = reverseBWT.length() - (reverseSuffixArray.get(alignment.loBound)+read.getReadLength()-alignment.gapOpens-alignment.gapExtensions) + 1;
                return Collections.<Alignment>singletonList(alignment);
            }

            //System.out.printf("Processing alignments; queue size = %d, alignment = %s, bound = %d%n", alignments.size(), alignment, lowerBounds.get(alignment.position+1).value);

            // if z < D(i) then return {}
            int mismatches = MAXIMUM_EDIT_DISTANCE - alignment.mismatches - alignment.gapOpens;            
            if( mismatches < lowerBounds.get(alignment.position+1).value )
                continue;

            if( alignment.mismatches > MAXIMUM_EDIT_DISTANCE )
                continue;

            if( alignment.state == AlignmentState.MATCH_MISMATCH ) {
                if( alignment.gapOpens < MAXIMUM_GAP_OPENS ) {
                    // Add a potential insertion.
                    BWAAlignment newAlignment = new BWAAlignment();
                    newAlignment.negativeStrand = alignment.negativeStrand;
                    newAlignment.position = alignment.position + 1;
                    newAlignment.state = AlignmentState.INSERTION;
                    newAlignment.gapOpens = alignment.gapOpens + 1;
                    newAlignment.gapExtensions = alignment.gapExtensions;

                    newAlignment.loBound = alignment.loBound;
                    newAlignment.hiBound = alignment.hiBound;

                    alignments.add(newAlignment);

                    // Add a potential deletion by marking a deletion and augmenting the position.
                    for(Base base: EnumSet.allOf(Base.class)) {
                        newAlignment = new BWAAlignment();
                        newAlignment.negativeStrand = alignment.negativeStrand;
                        newAlignment.position = alignment.position;
                        newAlignment.state = AlignmentState.DELETION;
                        newAlignment.gapOpens = alignment.gapOpens + 1;
                        newAlignment.gapExtensions = alignment.gapExtensions;

                        newAlignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
                        newAlignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);

                        if( newAlignment.loBound <= newAlignment.hiBound )
                            alignments.add(newAlignment);
                    }
                }
            }
            else if( alignment.state == AlignmentState.INSERTION ) {
                if( alignment.gapExtensions < MAXIMUM_GAP_EXTENSIONS ) {
                    // Add a potential insertion extension.
                    BWAAlignment newAlignment = new BWAAlignment();
                    newAlignment.negativeStrand = alignment.negativeStrand;
                    newAlignment.position = alignment.position + 1;
                    newAlignment.state = AlignmentState.INSERTION;
                    newAlignment.gapOpens = alignment.gapOpens;
                    newAlignment.gapExtensions = alignment.gapExtensions+1;

                    newAlignment.loBound = alignment.loBound;
                    newAlignment.hiBound = alignment.hiBound;                    

                    if( newAlignment.loBound <= newAlignment.hiBound )
                        alignments.add(newAlignment);
                }
            }
            else if( alignment.state == AlignmentState.DELETION ) {
                if( alignment.gapExtensions < MAXIMUM_GAP_EXTENSIONS ) {
                    // Add a potential deletion by marking a deletion and augmenting the position.
                    for(Base base: EnumSet.allOf(Base.class)) {
                        BWAAlignment newAlignment = new BWAAlignment();
                        newAlignment.negativeStrand = alignment.negativeStrand;
                        newAlignment.position = alignment.position;
                        newAlignment.state = AlignmentState.DELETION;
                        newAlignment.gapOpens = alignment.gapOpens;
                        newAlignment.gapExtensions = alignment.gapExtensions + 1;

                        newAlignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
                        newAlignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);

                        alignments.add(newAlignment);
                    }
                }
            }

            // Mismatches
            for(Base base: EnumSet.allOf(Base.class)) {
                // Create and initialize a new alignment, given that base as the candidate.
                BWAAlignment newAlignment = new BWAAlignment();
                newAlignment.negativeStrand = alignment.negativeStrand;
                newAlignment.position = alignment.position + 1;
                newAlignment.state = AlignmentState.MATCH_MISMATCH;
                newAlignment.mismatches = alignment.mismatches;
                if( base.toASCII() != bases[newAlignment.position] )
                    newAlignment.mismatches++;
                newAlignment.gapOpens = alignment.gapOpens;
                newAlignment.gapExtensions = alignment.gapExtensions;

                newAlignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
                newAlignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);

                // If this alignment is valid, add it to the list.
                if( newAlignment.loBound <= newAlignment.hiBound )
                    alignments.add(newAlignment);
            }
        }

        return Collections.emptyList();
    }

}
