package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.alignment.bwa.bwt.*;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.alignment.Aligner;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.File;
import java.util.*;

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
    private final int MAXIMUM_EDIT_DISTANCE = 4;

    /**
     * Maximum number of gap opens (-o option from original BWA).
     */
    private final int MAXIMUM_GAP_OPENS = 1;

    /**
     * Maximum number of gap extensions (-e option from original BWA).
     */
    private final int MAXIMUM_GAP_EXTENSIONS = 6;

    /**
     * Penalty for straight mismatches (-M option from original BWA).
     */
    public final int MISMATCH_PENALTY = 3;

    /**
     * Penalty for gap opens (-O option from original BWA).
     */
    public final int GAP_OPEN_PENALTY = 11;

    /**
     * Penalty for gap extensions (-E option from original BWA).
     */
    public final int GAP_EXTENSION_PENALTY = 4;

    public BWAAligner( File forwardBWTFile, File reverseBWTFile, File reverseSuffixArrayFile ) {
        forwardBWT = new BWTReader(forwardBWTFile).read();
        reverseBWT = new BWTReader(reverseBWTFile).read();
        reverseSuffixArray = new SuffixArrayReader(reverseSuffixArrayFile).read();
    }

    public List<Alignment> align( SAMRecord read ) {
        byte[] forwardBases = read.getReadBases();
        byte[] reverseBases = BaseUtils.simpleReverseComplement(forwardBases);

        List<LowerBound> forwardLowerBounds = LowerBound.create(forwardBases,forwardBWT);
        List<LowerBound> reverseLowerBounds = LowerBound.create(reverseBases,reverseBWT);

        //for( int i = 0; i < forwardLowerBounds.size(); i++ )
        //    System.out.printf("ForwardBWT: lb[%d] = %s%n",i,forwardLowerBounds.get(i));
        //for( int i = 0; i < reverseLowerBounds.size(); i++ )
        //    System.out.printf("ReverseBWT: lb[%d] = %s%n",i,reverseLowerBounds.get(i));

        PriorityQueue<BWAAlignment> alignments = new PriorityQueue<BWAAlignment>();

        // Create a fictional initial alignment, with the position just off the end of the read, and the limits
        // set as the entire BWT.
        alignments.add(createSeedAlignment(forwardBWT));
        //alignments.add(createSeedAlignment(reverseBWT));

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
                    // Add a potential insertion extension.
                    BWAAlignment insertionAlignment = createInsertionAlignment(alignment);
                    insertionAlignment.gapOpens++;
                    alignments.add(insertionAlignment);

                    // Add a potential deletion by marking a deletion and augmenting the position.
                    List<BWAAlignment> newAlignments = createDeletionAlignments(bwt,alignment);
                    for( BWAAlignment deletionAlignment: newAlignments )
                        deletionAlignment.gapOpens++;
                    alignments.addAll(newAlignments);
                }
            }
            else if( alignment.state == AlignmentState.INSERTION ) {
                if( alignment.gapExtensions < MAXIMUM_GAP_EXTENSIONS ) {
                    // Add a potential insertion extension.
                    BWAAlignment newAlignment = createInsertionAlignment(alignment);
                    newAlignment.gapExtensions++;
                    alignments.add(newAlignment);
                }
            }
            else if( alignment.state == AlignmentState.DELETION ) {
                if( alignment.gapExtensions < MAXIMUM_GAP_EXTENSIONS ) {
                    // Add a potential deletion by marking a deletion and augmenting the position.
                    List<BWAAlignment> newAlignments = createDeletionAlignments(bwt,alignment);
                    for( BWAAlignment newAlignment: newAlignments )
                        newAlignment.gapExtensions++;
                    alignments.addAll(newAlignments);
                }
            }

            // Mismatches
            alignments.addAll(createMatchedAlignments(bwt,alignment,bases));
        }

        return Collections.emptyList();
    }

    /**
     * Create an seeding alignment to use as a starting point when traversing.
     * @param bwt source BWT.
     * @return Seed alignment.
     */
    private BWAAlignment createSeedAlignment(BWT bwt) {
        BWAAlignment seed = new BWAAlignment(this);
        seed.negativeStrand = (bwt == forwardBWT);
        seed.position = 0;
        seed.loBound = 0;
        seed.hiBound = bwt.length();
        seed.state = AlignmentState.MATCH_MISMATCH;
        seed.mismatches = 0;
        return seed;
    }

    /**
     * Creates a new alignments representing direct matches / mismatches.
     * @param bwt Source BWT with which to work.
     * @param alignment Alignment for the previous position.
     * @param bases The bases in the read.
     * @return New alignment representing this position if valid; null otherwise.
     */
    private List<BWAAlignment> createMatchedAlignments( BWT bwt, BWAAlignment alignment, byte[] bases ) {
        List<BWAAlignment> newAlignments = new ArrayList<BWAAlignment>();
        for(Base base: EnumSet.allOf(Base.class)) {
            BWAAlignment newAlignment = alignment.clone();

            newAlignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
            newAlignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);

            // If this alignment is valid, skip it.
            if( newAlignment.loBound > newAlignment.hiBound )
                continue;

            newAlignment.position++;
            newAlignment.state = AlignmentState.MATCH_MISMATCH;
            if( base.toASCII() != bases[newAlignment.position] )
                newAlignment.mismatches++;

            newAlignments.add(newAlignment);
        }

        return newAlignments;
    }

    /**
     * Create a new alignment representing an insertion at this point in the read.
     * @param alignment Alignment from which to derive the insertion.
     * @return New alignment reflecting the insertion.
     */
    private BWAAlignment createInsertionAlignment( BWAAlignment alignment ) {
        // Add a potential insertion extension.
        BWAAlignment newAlignment = alignment.clone();
        newAlignment.position++;
        newAlignment.state = AlignmentState.INSERTION;
        return newAlignment;
    }

    /**
     * Create new alignments representing a deletion a this point in the read.
     * @param alignment Alignment from which to derive the deletion.
     * @return New alignments reflecting all possible deletions.
     */
    private List<BWAAlignment> createDeletionAlignments( BWT bwt, BWAAlignment alignment) {
        List<BWAAlignment> newAlignments = new ArrayList<BWAAlignment>();
        for(Base base: EnumSet.allOf(Base.class)) {
            BWAAlignment newAlignment = alignment.clone();

            newAlignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
            newAlignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);

            // If this alignment is valid, skip it.
            if( newAlignment.loBound > newAlignment.hiBound )
                continue;

            newAlignment.state = AlignmentState.DELETION;            

            newAlignments.add(newAlignment);
        }

        return newAlignments;
    }

}
