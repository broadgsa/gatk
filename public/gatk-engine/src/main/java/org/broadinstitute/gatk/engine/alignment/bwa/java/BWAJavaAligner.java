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

package org.broadinstitute.gatk.engine.alignment.bwa.java;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.engine.alignment.Alignment;
import org.broadinstitute.gatk.engine.alignment.bwa.BWAAligner;
import org.broadinstitute.gatk.engine.alignment.bwa.BWAConfiguration;
import org.broadinstitute.gatk.engine.alignment.reference.bwt.*;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Create imperfect alignments from the read to the genome represented by the given BWT / suffix array. 
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAJavaAligner extends BWAAligner {
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
    private SuffixArray forwardSuffixArray;

    /**
     * Suffix array in the reverse direction.
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

    /**
     * Skip the ends of indels.
     */
    public final int INDEL_END_SKIP = 5;

    public BWAJavaAligner( File forwardBWTFile, File reverseBWTFile, File forwardSuffixArrayFile, File reverseSuffixArrayFile ) {
        super(null,null);
        forwardBWT = new BWTReader(forwardBWTFile).read();
        reverseBWT = new BWTReader(reverseBWTFile).read();
        forwardSuffixArray = new SuffixArrayReader(forwardSuffixArrayFile,forwardBWT).read();
        reverseSuffixArray = new SuffixArrayReader(reverseSuffixArrayFile,reverseBWT).read();
    }

    /**
     * Close this instance of the BWA pointer and delete its resources.
     */
    @Override
    public void close()  {
        throw new UnsupportedOperationException("BWA aligner can't currently be closed.");
    }

    /**
     * Update the current parameters of this aligner.
     * @param configuration New configuration to set.
     */
    public void updateConfiguration(BWAConfiguration configuration) {
        throw new UnsupportedOperationException("Configuration of the BWA aligner can't currently be changed.");
    }

    /**
     * Allow the aligner to choose one alignment randomly from the pile of best alignments.
     * @param bases Bases to align.
     * @return An align
     */
    public Alignment getBestAlignment(final byte[] bases) { throw new UnsupportedOperationException("BWAJavaAligner does not yet support the standard Aligner interface."); }

    /**
     * Align the read to the reference.
     * @param read Read to align.
     * @param header Optional header to drop in place.
     * @return A list of the alignments.
     */
    public SAMRecord align(final SAMRecord read, final SAMFileHeader header) { throw new UnsupportedOperationException("BWAJavaAligner does not yet support the standard Aligner interface."); }

    /**
     * Get a iterator of alignments, batched by mapping quality.
     * @param bases List of bases.
     * @return Iterator to alignments.
     */
    public Iterable<Alignment[]> getAllAlignments(final byte[] bases) { throw new UnsupportedOperationException("BWAJavaAligner does not yet support the standard Aligner interface."); }

    /**
     * Get a iterator of aligned reads, batched by mapping quality.
     * @param read Read to align.
     * @param newHeader Optional new header to use when aligning the read.  If present, it must be null.
     * @return Iterator to alignments.
     */
    public Iterable<SAMRecord[]> alignAll(final SAMRecord read, final SAMFileHeader newHeader) { throw new UnsupportedOperationException("BWAJavaAligner does not yet support the standard Aligner interface."); }


    public List<Alignment> align( SAMRecord read ) {
        List<Alignment> successfulMatches = new ArrayList<Alignment>();

        Byte[] uncomplementedBases = normalizeBases(read.getReadBases());
        Byte[] complementedBases = normalizeBases(Utils.reverse(BaseUtils.simpleReverseComplement(read.getReadBases())));

        List<LowerBound> forwardLowerBounds = LowerBound.create(uncomplementedBases,forwardBWT);
        List<LowerBound> reverseLowerBounds = LowerBound.create(complementedBases,reverseBWT);

        // Seed the best score with any score that won't overflow on comparison.
        int bestScore = Integer.MAX_VALUE - MISMATCH_PENALTY;
        int bestDiff = MAXIMUM_EDIT_DISTANCE+1;
        int maxDiff = MAXIMUM_EDIT_DISTANCE;

        PriorityQueue<BWAAlignment> alignments = new PriorityQueue<BWAAlignment>();

        // Create a fictional initial alignment, with the position just off the end of the read, and the limits
        // set as the entire BWT.
        alignments.add(createSeedAlignment(reverseBWT));
        alignments.add(createSeedAlignment(forwardBWT));

        while(!alignments.isEmpty()) {
            BWAAlignment alignment = alignments.remove();

            // From bwtgap.c in the original BWT; if the rank is worse than the best score + the mismatch PENALTY, move on.
            if( alignment.getScore() > bestScore + MISMATCH_PENALTY )
                break;

            Byte[] bases = alignment.isNegativeStrand() ? complementedBases : uncomplementedBases;
            BWT bwt = alignment.isNegativeStrand() ? forwardBWT : reverseBWT;
            List<LowerBound> lowerBounds = alignment.isNegativeStrand() ? reverseLowerBounds : forwardLowerBounds;

            // if z < D(i) then return {}
            int mismatches = maxDiff - alignment.getMismatches() - alignment.getGapOpens() - alignment.getGapExtensions();
            if( alignment.position < lowerBounds.size()-1 && mismatches < lowerBounds.get(alignment.position+1).value )
                continue;

            if(mismatches == 0) {
                exactMatch(alignment,bases,bwt);
                if(alignment.loBound > alignment.hiBound)
                    continue;
            }

            // Found a valid alignment; store it and move on.
            if(alignment.position >= read.getReadLength()-1) {
                for(long bwtIndex = alignment.loBound; bwtIndex <= alignment.hiBound; bwtIndex++) {
                    BWAAlignment finalAlignment = alignment.clone();

                    if( finalAlignment.isNegativeStrand() )
                        finalAlignment.setAlignmentStart(forwardSuffixArray.get(bwtIndex) + 1);
                    else {
                        int sizeAlongReference = read.getReadLength() -
                                finalAlignment.getNumberOfBasesMatchingState(AlignmentState.INSERTION) +
                                finalAlignment.getNumberOfBasesMatchingState(AlignmentState.DELETION);
                        finalAlignment.setAlignmentStart(reverseBWT.length() - reverseSuffixArray.get(bwtIndex) - sizeAlongReference + 1);
                    }

                    successfulMatches.add(finalAlignment);

                    bestScore = Math.min(finalAlignment.getScore(),bestScore);
                    bestDiff = Math.min(finalAlignment.getMismatches()+finalAlignment.getGapOpens()+finalAlignment.getGapExtensions(),bestDiff);
                    maxDiff = bestDiff + 1;
                }

                continue;
            }

            //System.out.printf("Processing alignments; queue size = %d, alignment = %s, bound = %d, base = %s%n", alignments.size(), alignment, lowerBounds.get(alignment.position+1).value, alignment.position >= 0 ? (char)bases[alignment.position].byteValue() : "");
            /*
            System.out.printf("#1\t[%d,%d,%d,%c]\t[%d,%d,%d]\t[%d,%d]\t[%d,%d]%n",alignments.size(),
                                                        alignment.negativeStrand?1:0,
                                                        bases.length-alignment.position-1,
                                                        alignment.getCurrentState().toString().charAt(0),
                                                        alignment.getMismatches(),
                                                        alignment.getGapOpens(),
                                                        alignment.getGapExtensions(),
                                                        lowerBounds.get(alignment.position+1).value,
                                                        lowerBounds.get(alignment.position+1).width,
                                                        alignment.loBound,
                                                        alignment.hiBound);
                                                        */

            // Temporary -- look ahead to see if the next alignment is bounded.
            boolean allowDifferences = mismatches > 0;
            boolean allowMismatches = mismatches > 0;

            if( allowDifferences &&
                alignment.position+1 >= INDEL_END_SKIP-1+alignment.getGapOpens()+alignment.getGapExtensions() &&
                read.getReadLength()-1-(alignment.position+1) >= INDEL_END_SKIP+alignment.getGapOpens()+alignment.getGapExtensions() ) {
                if( alignment.getCurrentState() == AlignmentState.MATCH_MISMATCH ) {
                    if( alignment.getGapOpens() < MAXIMUM_GAP_OPENS ) {
                        // Add a potential insertion extension.
                        BWAAlignment insertionAlignment = createInsertionAlignment(alignment);
                        insertionAlignment.incrementGapOpens();
                        alignments.add(insertionAlignment);

                        // Add a potential deletion by marking a deletion and augmenting the position.
                        List<BWAAlignment> deletionAlignments = createDeletionAlignments(bwt,alignment);
                        for( BWAAlignment deletionAlignment: deletionAlignments )
                            deletionAlignment.incrementGapOpens();
                        alignments.addAll(deletionAlignments);
                    }
                }
                else if( alignment.getCurrentState() == AlignmentState.INSERTION ) {
                    if( alignment.getGapExtensions() < MAXIMUM_GAP_EXTENSIONS && mismatches > 0 ) {
                        // Add a potential insertion extension.
                        BWAAlignment insertionAlignment = createInsertionAlignment(alignment);
                        insertionAlignment.incrementGapExtensions();
                        alignments.add(insertionAlignment);
                    }
                }
                else if( alignment.getCurrentState() == AlignmentState.DELETION ) {
                    if( alignment.getGapExtensions() < MAXIMUM_GAP_EXTENSIONS && mismatches > 0 ) {
                        // Add a potential deletion by marking a deletion and augmenting the position.
                        List<BWAAlignment> deletionAlignments = createDeletionAlignments(bwt,alignment);
                        for( BWAAlignment deletionAlignment: deletionAlignments )
                            deletionAlignment.incrementGapExtensions();
                        alignments.addAll(deletionAlignments);
                    }
                }
            }

            // Mismatches
            alignments.addAll(createMatchedAlignments(bwt,alignment,bases,allowDifferences&&allowMismatches));
        }

        return successfulMatches;
    }

    /**
     * Create an seeding alignment to use as a starting point when traversing.
     * @param bwt source BWT.
     * @return Seed alignment.
     */
    private BWAAlignment createSeedAlignment(BWT bwt) {
        BWAAlignment seed = new BWAAlignment(this);
        seed.setNegativeStrand(bwt == forwardBWT);
        seed.position = -1;
        seed.loBound = 0;
        seed.hiBound = bwt.length();
        return seed;
    }

    /**
     * Creates a new alignments representing direct matches / mismatches.
     * @param bwt Source BWT with which to work.
     * @param alignment Alignment for the previous position.
     * @param bases The bases in the read.
     * @param allowMismatch Should mismatching bases be allowed?
     * @return New alignment representing this position if valid; null otherwise.
     */
    private List<BWAAlignment> createMatchedAlignments( BWT bwt, BWAAlignment alignment, Byte[] bases, boolean allowMismatch ) {
        List<BWAAlignment> newAlignments = new ArrayList<BWAAlignment>();

        List<Byte> baseChoices = new ArrayList<Byte>();
        Byte thisBase = bases[alignment.position+1];

        if( allowMismatch )
            baseChoices.addAll(Bases.allOf());
        else
            baseChoices.add(thisBase);

        if( thisBase != null ) {
            // Keep rotating the current base to the last position until we've hit the current base.
            for( ;; ) {
                baseChoices.add(baseChoices.remove(0));
                if( thisBase.equals(baseChoices.get(baseChoices.size()-1)) )
                    break;

            }
        }

        for(byte base: baseChoices) {
            BWAAlignment newAlignment = alignment.clone();

            newAlignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
            newAlignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);

            // If this alignment is valid, skip it.
            if( newAlignment.loBound > newAlignment.hiBound )
                continue;

            newAlignment.position++;
            newAlignment.addState(AlignmentState.MATCH_MISMATCH);
            if( bases[newAlignment.position] == null || base != bases[newAlignment.position] )
                newAlignment.incrementMismatches();

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
        newAlignment.addState(AlignmentState.INSERTION);
        return newAlignment;
    }

    /**
     * Create new alignments representing a deletion at this point in the read.
     * @param bwt source BWT for inferring deletion info.
     * @param alignment Alignment from which to derive the deletion.
     * @return New alignments reflecting all possible deletions.
     */
    private List<BWAAlignment> createDeletionAlignments( BWT bwt, BWAAlignment alignment) {
        List<BWAAlignment> newAlignments = new ArrayList<BWAAlignment>();
        for(byte base: Bases.instance) {
            BWAAlignment newAlignment = alignment.clone();

            newAlignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
            newAlignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);

            // If this alignment is valid, skip it.
            if( newAlignment.loBound > newAlignment.hiBound )
                continue;

            newAlignment.addState(AlignmentState.DELETION);

            newAlignments.add(newAlignment);
        }

        return newAlignments;
    }

    /**
     * Exactly match the given alignment against the given BWT.
     * @param alignment Alignment to match.
     * @param bases Bases to use.
     * @param bwt BWT to use.
     */
    private void exactMatch( BWAAlignment alignment, Byte[] bases, BWT bwt ) {
        while( ++alignment.position < bases.length ) {
            byte base = bases[alignment.position];
            alignment.loBound = bwt.counts(base) + bwt.occurrences(base,alignment.loBound-1) + 1;
            alignment.hiBound = bwt.counts(base) + bwt.occurrences(base,alignment.hiBound);
            if( alignment.loBound > alignment.hiBound )
                return;
        }
    }

    /**
     * Make each base into A/C/G/T or null if unknown.
     * @param bases Base string to normalize.
     * @return Array of normalized bases.
     */
    private Byte[] normalizeBases( byte[] bases ) {
        Byte[] normalBases = new Byte[bases.length];
        for(int i = 0; i < bases.length; i++)
            normalBases[i] = Bases.fromASCII(bases[i]);
        return normalBases;
    }
}
