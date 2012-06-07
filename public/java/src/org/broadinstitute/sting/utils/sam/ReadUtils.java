/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.utils.sam;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.NGSPlatform;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.util.*;

/**
 * A miscellaneous collection of utilities for working with SAM files, headers, etc.
 * Static methods only, please.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadUtils {
    
    private static final String OFFSET_OUT_OF_BOUNDS_EXCEPTION = "Offset cannot be greater than read length %d : %d";
    private static final String OFFSET_NOT_ZERO_EXCEPTION = "We ran past the end of the read and never found the offset, something went wrong!";
    
    private ReadUtils() {
    }

    private static int DEFAULT_ADAPTOR_SIZE = 100;
    public static int CLIPPING_GOAL_NOT_REACHED = -1;

    /**
     * A marker to tell which end of the read has been clipped
     */
    public enum ClippingTail {
        LEFT_TAIL,
        RIGHT_TAIL
    }

    /**
     * A HashMap of the SAM spec read flag names
     *
     * Note: This is not being used right now, but can be useful in the future
     */
    private static final Map<Integer, String> readFlagNames = new HashMap<Integer, String>();

    static {
        readFlagNames.put(0x1, "Paired");
        readFlagNames.put(0x2, "Proper");
        readFlagNames.put(0x4, "Unmapped");
        readFlagNames.put(0x8, "MateUnmapped");
        readFlagNames.put(0x10, "Forward");
        //readFlagNames.put(0x20, "MateForward");
        readFlagNames.put(0x40, "FirstOfPair");
        readFlagNames.put(0x80, "SecondOfPair");
        readFlagNames.put(0x100, "NotPrimary");
        readFlagNames.put(0x200, "NON-PF");
        readFlagNames.put(0x400, "Duplicate");
    }

    /**
     * This enum represents all the different ways in which a read can overlap an interval.
     *
     * NO_OVERLAP_CONTIG:
     * read and interval are in different contigs.
     *
     * NO_OVERLAP_LEFT:
     * the read does not overlap the interval.
     *
     *                        |----------------| (interval)
     *   <---------------->                      (read)
     *
     * NO_OVERLAP_RIGHT:
     * the read does not overlap the interval.
     *
     *   |----------------|                      (interval)
     *                        <----------------> (read)
     *
     * OVERLAP_LEFT:
     * the read starts before the beginning of the interval but ends inside of it
     *
     *          |----------------| (interval)
     *   <---------------->        (read)
     *
     * OVERLAP_RIGHT:
     * the read starts inside the interval but ends outside of it
     *
     *   |----------------|     (interval)
     *       <----------------> (read)
     *
     * OVERLAP_LEFT_AND_RIGHT:
     * the read starts before the interval and ends after the interval
     *
     *      |-----------|     (interval)
     *  <-------------------> (read)
     *
     * OVERLAP_CONTAINED:
     * the read starts and ends inside the interval
     *
     *  |----------------|     (interval)
     *     <-------->          (read)
     */
    public enum ReadAndIntervalOverlap {NO_OVERLAP_CONTIG, NO_OVERLAP_LEFT, NO_OVERLAP_RIGHT, NO_OVERLAP_HARDCLIPPED_LEFT, NO_OVERLAP_HARDCLIPPED_RIGHT, OVERLAP_LEFT, OVERLAP_RIGHT, OVERLAP_LEFT_AND_RIGHT, OVERLAP_CONTAINED}

    /**
     * Creates a SAMFileWriter with the given compression level if you request a bam file. Creates a regular
     * SAMFileWriter without compression otherwise.
     *
     * @param header
     * @param presorted
     * @param file
     * @param compression
     * @return a SAMFileWriter with the compression level if it is a bam.
     */
    public static SAMFileWriter createSAMFileWriterWithCompression(SAMFileHeader header, boolean presorted, String file, int compression) {
        if (file.endsWith(".bam"))
            return new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(file), compression);
        return new SAMFileWriterFactory().makeSAMOrBAMWriter(header, presorted, new File(file));
    }

    /**
     * is this base inside the adaptor of the read?
     *
     * There are two cases to treat here:
     *
     * 1) Read is in the negative strand => Adaptor boundary is on the left tail
     * 2) Read is in the positive strand => Adaptor boundary is on the right tail
     *
     * Note: We return false to all reads that are UNMAPPED or have an weird big insert size (probably due to mismapping or bigger event)
     *
     * @param read the read to test
     * @param basePos base position in REFERENCE coordinates (not read coordinates)
     * @return whether or not the base is in the adaptor
     */
    public static boolean isBaseInsideAdaptor(final GATKSAMRecord read, long basePos) {
        Integer adaptorBoundary = getAdaptorBoundary(read);
        if (adaptorBoundary == null || read.getInferredInsertSize() > DEFAULT_ADAPTOR_SIZE)
            return false;

        return read.getReadNegativeStrandFlag() ? basePos <= adaptorBoundary : basePos >= adaptorBoundary;
    }

    /**
     * Finds the adaptor boundary around the read and returns the first base inside the adaptor that is closest to
     * the read boundary. If the read is in the positive strand, this is the first base after the end of the
     * fragment (Picard calls it 'insert'), if the read is in the negative strand, this is the first base before the
     * beginning of the fragment.
     *
     * There are two cases we need to treat here:
     *
     * 1) Our read is in the reverse strand :
     *
     *     <----------------------| *
     *   |--------------------->
     *
     *   in these cases, the adaptor boundary is at the mate start (minus one)
     *
     * 2) Our read is in the forward strand :
     *
     *   |---------------------->   *
     *     <----------------------|
     *
     *   in these cases the adaptor boundary is at the start of the read plus the inferred insert size (plus one)
     *
     * @param read the read being tested for the adaptor boundary
     * @return the reference coordinate for the adaptor boundary (effectively the first base IN the adaptor, closest to the read. NULL if the read is unmapped or the mate is mapped to another contig.
     */
    public static Integer getAdaptorBoundary(final SAMRecord read) {
        final int MAXIMUM_ADAPTOR_LENGTH = 8;
        final int insertSize = Math.abs(read.getInferredInsertSize());    // the inferred insert size can be negative if the mate is mapped before the read (so we take the absolute value)

        if (insertSize == 0 || read.getReadUnmappedFlag())                // no adaptors in reads with mates in another chromosome or unmapped pairs
            return null;                                                  
        
        Integer adaptorBoundary;                                          // the reference coordinate for the adaptor boundary (effectively the first base IN the adaptor, closest to the read)
        if (read.getReadNegativeStrandFlag())
            adaptorBoundary = read.getMateAlignmentStart() - 1;           // case 1 (see header)
        else
            adaptorBoundary = read.getAlignmentStart() + insertSize + 1;  // case 2 (see header)

        if ( (adaptorBoundary < read.getAlignmentStart() - MAXIMUM_ADAPTOR_LENGTH) || (adaptorBoundary > read.getAlignmentEnd() + MAXIMUM_ADAPTOR_LENGTH) )
            adaptorBoundary = null;                                       // we are being conservative by not allowing the adaptor boundary to go beyond what we belive is the maximum size of an adaptor
        
        return adaptorBoundary;
    }

    /**
     * is the read a 454 read?
     *
     * @param read the read to test
     * @return checks the read group tag PL for the default 454 tag
     */
    public static boolean is454Read(SAMRecord read) {
        return NGSPlatform.fromRead(read) == NGSPlatform.LS454;
    }

    /**
     * is the read an IonTorrent read?
     *
     * @param read the read to test
     * @return checks the read group tag PL for the default ion tag
     */
    public static boolean isIonRead(SAMRecord read) {
        return NGSPlatform.fromRead(read) == NGSPlatform.ION_TORRENT;
    }

    /**
     * is the read a SOLiD read?
     *
     * @param read the read to test
     * @return checks the read group tag PL for the default SOLiD tag
     */
    public static boolean isSOLiDRead(SAMRecord read) {
        return NGSPlatform.fromRead(read) == NGSPlatform.SOLID;
    }

    /**
     * is the read a SLX read?
     *
     * @param read the read to test
     * @return checks the read group tag PL for the default SLX tag
     */
    public static boolean isIlluminaRead(SAMRecord read) {
        return NGSPlatform.fromRead(read) == NGSPlatform.ILLUMINA;
    }

    /**
     * checks if the read has a platform tag in the readgroup equal to 'name'.
     * Assumes that 'name' is upper-cased.
     *
     * @param read the read to test
     * @param name the upper-cased platform name to test
     * @return whether or not name == PL tag in the read group of read
     */
    public static boolean isPlatformRead(SAMRecord read, String name) {

        SAMReadGroupRecord readGroup = read.getReadGroup();
        if (readGroup != null) {
            Object readPlatformAttr = readGroup.getAttribute("PL");
            if (readPlatformAttr != null)
                return readPlatformAttr.toString().toUpperCase().contains(name);
        }
        return false;
    }


    /**
     * Returns the collections of reads sorted in coordinate order, according to the order defined
     * in the reads themselves
     *
     * @param reads
     * @return
     */
    public final static List<GATKSAMRecord> sortReadsByCoordinate(List<GATKSAMRecord> reads) {
        final SAMRecordComparator comparer = new SAMRecordCoordinateComparator();
        Collections.sort(reads, comparer);
        return reads;
    }

    /**
     * If a read starts in INSERTION, returns the first element length.
     *
     * Warning: If the read has Hard or Soft clips before the insertion this function will return 0.
     *
     * @param read
     * @return the length of the first insertion, or 0 if there is none (see warning).
     */
    public final static int getFirstInsertionOffset(SAMRecord read) {
        CigarElement e = read.getCigar().getCigarElement(0);
        if ( e.getOperator() == CigarOperator.I )
            return e.getLength();
        else
            return 0;
    }

    /**
     * If a read ends in INSERTION, returns the last element length.
     *
     * Warning: If the read has Hard or Soft clips after the insertion this function will return 0.
     *
     * @param read
     * @return the length of the last insertion, or 0 if there is none (see warning).
     */
    public final static int getLastInsertionOffset(SAMRecord read) {
        CigarElement e = read.getCigar().getCigarElement(read.getCigarLength() - 1);
        if ( e.getOperator() == CigarOperator.I )
            return e.getLength();
        else
            return 0;
    }

    /**
     * Determines what is the position of the read in relation to the interval.
     * Note: This function uses the UNCLIPPED ENDS of the reads for the comparison.
     * @param read the read
     * @param interval the interval
     * @return the overlap type as described by ReadAndIntervalOverlap enum (see above)
     */
    public static ReadAndIntervalOverlap getReadAndIntervalOverlapType(GATKSAMRecord read, GenomeLoc interval) {

        int sStart = read.getSoftStart();
        int sStop = read.getSoftEnd();
        int uStart = read.getUnclippedStart();
        int uStop = read.getUnclippedEnd();

        if ( !read.getReferenceName().equals(interval.getContig()) )
            return ReadAndIntervalOverlap.NO_OVERLAP_CONTIG;

        else if ( uStop < interval.getStart() )
            return ReadAndIntervalOverlap.NO_OVERLAP_LEFT;

        else if ( uStart > interval.getStop() )
            return ReadAndIntervalOverlap.NO_OVERLAP_RIGHT;

        else if ( sStop < interval.getStart() )
            return ReadAndIntervalOverlap.NO_OVERLAP_HARDCLIPPED_LEFT;

        else if ( sStart > interval.getStop() )
            return ReadAndIntervalOverlap.NO_OVERLAP_HARDCLIPPED_RIGHT;

        else if ( (sStart >= interval.getStart()) &&
                  (sStop <= interval.getStop()) )
            return ReadAndIntervalOverlap.OVERLAP_CONTAINED;

        else if ( (sStart < interval.getStart()) &&
                  (sStop > interval.getStop()) )
            return ReadAndIntervalOverlap.OVERLAP_LEFT_AND_RIGHT;

        else if ( (sStart < interval.getStart()) )
            return ReadAndIntervalOverlap.OVERLAP_LEFT;

        else
            return ReadAndIntervalOverlap.OVERLAP_RIGHT;
    }

    /**
     * Pre-processes the results of getReadCoordinateForReferenceCoordinate(GATKSAMRecord, int) to take care of
     * two corner cases:
     * 
     * 1. If clipping the right tail (end of the read) getReadCoordinateForReferenceCoordinate and fall inside
     * a deletion return the base after the deletion. If clipping the left tail (beginning of the read) it
     * doesn't matter because it already returns the previous base by default.
     * 
     * 2. If clipping the left tail (beginning of the read) getReadCoordinateForReferenceCoordinate and the
     * read starts with an insertion, and you're requesting the first read based coordinate, it will skip
     * the leading insertion (because it has the same reference coordinate as the following base).
     *
     * @param read
     * @param refCoord
     * @param tail
     * @return the read coordinate corresponding to the requested reference coordinate for clipping.
     */
    @Requires({"refCoord >= read.getUnclippedStart()", "refCoord <= read.getUnclippedEnd() || (read.getUnclippedEnd() < read.getUnclippedStart())"})
    @Ensures({"result >= 0", "result < read.getReadLength()"})
    public static int getReadCoordinateForReferenceCoordinate(GATKSAMRecord read, int refCoord, ClippingTail tail) {
        return getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), refCoord, tail, false);
    }

    public static int getReadCoordinateForReferenceCoordinate(final int alignmentStart, final Cigar cigar, final int refCoord, final ClippingTail tail, final boolean allowGoalNotReached) {
        Pair<Integer, Boolean> result = getReadCoordinateForReferenceCoordinate(alignmentStart, cigar, refCoord, allowGoalNotReached);
        int readCoord = result.getFirst();

        // Corner case one: clipping the right tail and falls on deletion, move to the next
        // read coordinate. It is not a problem for the left tail because the default answer
        // from getReadCoordinateForReferenceCoordinate is to give the previous read coordinate.
        if (result.getSecond() && tail == ClippingTail.RIGHT_TAIL)
            readCoord++;

        // clipping the left tail and first base is insertion, go to the next read coordinate
        // with the same reference coordinate. Advance to the next cigar element, or to the
        // end of the read if there is no next element.
        Pair<Boolean, CigarElement> firstElementIsInsertion = readStartsWithInsertion(cigar);
        if (readCoord == 0 && tail == ClippingTail.LEFT_TAIL && firstElementIsInsertion.getFirst())
            readCoord = Math.min(firstElementIsInsertion.getSecond().getLength(), cigar.getReadLength() - 1);

        return readCoord;
    }

    /**
     * Returns the read coordinate corresponding to the requested reference coordinate.
     *
     * WARNING: if the requested reference coordinate happens to fall inside a deletion in the read, this function
     * will return the last read base before the deletion. This function returns a
     * Pair(int readCoord, boolean fallsInsideDeletion) so you can choose which readCoordinate to use when faced with
     * a deletion.
     *
     * SUGGESTION: Use getReadCoordinateForReferenceCoordinate(GATKSAMRecord, int, ClippingTail) instead to get a
     * pre-processed result according to normal clipping needs. Or you can use this function and tailor the
     * behavior to your needs.
     *
     * @param read
     * @param refCoord
     * @return the read coordinate corresponding to the requested reference coordinate. (see warning!)
     */
    @Requires({"refCoord >= read.getSoftStart()", "refCoord <= read.getSoftEnd()"})
    @Ensures({"result.getFirst() >= 0", "result.getFirst() < read.getReadLength()"})
    public static Pair<Integer, Boolean> getReadCoordinateForReferenceCoordinate(GATKSAMRecord read, int refCoord) {
        return getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), refCoord, false);
    }

    public static Pair<Integer, Boolean> getReadCoordinateForReferenceCoordinate(final int alignmentStart, final Cigar cigar, final int refCoord, final boolean allowGoalNotReached) {
        int readBases = 0;
        int refBases = 0;
        boolean fallsInsideDeletion = false;

        int goal = refCoord - alignmentStart;  // The goal is to move this many reference bases
        if (goal < 0) {
            if (allowGoalNotReached) {
                return new Pair<Integer, Boolean>(CLIPPING_GOAL_NOT_REACHED, false);
            } else {
                throw new ReviewedStingException("Somehow the requested coordinate is not covered by the read. Too many deletions?");
            }
        }
        boolean goalReached = refBases == goal;

        Iterator<CigarElement> cigarElementIterator = cigar.getCigarElements().iterator();
        while (!goalReached && cigarElementIterator.hasNext()) {
            CigarElement cigarElement = cigarElementIterator.next();
            int shift = 0;

            if (cigarElement.getOperator().consumesReferenceBases() || cigarElement.getOperator() == CigarOperator.SOFT_CLIP) {
                if (refBases + cigarElement.getLength() < goal)
                    shift = cigarElement.getLength();
                else
                    shift = goal - refBases;

                refBases += shift;
            }
            goalReached = refBases == goal;

            if (!goalReached && cigarElement.getOperator().consumesReadBases())
                readBases += cigarElement.getLength();

            if (goalReached) {
                // Is this base's reference position within this cigar element? Or did we use it all?
                boolean endsWithinCigar = shift < cigarElement.getLength();

                // If it isn't, we need to check the next one. There should *ALWAYS* be a next one
                // since we checked if the goal coordinate is within the read length, so this is just a sanity check.
                if (!endsWithinCigar && !cigarElementIterator.hasNext()) {
                    if (allowGoalNotReached) {
                        return new Pair<Integer, Boolean>(CLIPPING_GOAL_NOT_REACHED, false);
                    } else {
                        throw new ReviewedStingException("Reference coordinate corresponds to a non-existent base in the read. This should never happen -- call Mauricio");
                    }
                }

                CigarElement nextCigarElement;

                // if we end inside the current cigar element, we just have to check if it is a deletion
                if (endsWithinCigar)
                    fallsInsideDeletion = cigarElement.getOperator() == CigarOperator.DELETION;

                // if we end outside the current cigar element, we need to check if the next element is an insertion or deletion.
                else {
                    nextCigarElement = cigarElementIterator.next();

                    // if it's an insertion, we need to clip the whole insertion before looking at the next element
                    if (nextCigarElement.getOperator() == CigarOperator.INSERTION) {
                        readBases += nextCigarElement.getLength();
                        if (!cigarElementIterator.hasNext()) {
                            if (allowGoalNotReached) {
                                return new Pair<Integer, Boolean>(CLIPPING_GOAL_NOT_REACHED, false);
                            } else {
                                throw new ReviewedStingException("Reference coordinate corresponds to a non-existent base in the read. This should never happen -- call Mauricio");
                            }
                        }

                        nextCigarElement = cigarElementIterator.next();
                    }

                    // if it's a deletion, we will pass the information on to be handled downstream.
                    fallsInsideDeletion = nextCigarElement.getOperator() == CigarOperator.DELETION;
                }

                // If we reached our goal outside a deletion, add the shift
                if (!fallsInsideDeletion && cigarElement.getOperator().consumesReadBases())
                    readBases += shift;

                // If we reached our goal inside a deletion, but the deletion is the next cigar element then we need
                // to add the shift of the current cigar element but go back to it's last element to return the last
                // base before the deletion (see warning in function contracts)
                else if (fallsInsideDeletion && !endsWithinCigar)
                    readBases += shift - 1;

                // If we reached our goal inside a deletion then we must backtrack to the last base before the deletion
                else if (fallsInsideDeletion && endsWithinCigar)
                    readBases--;
            }
        }

        if (!goalReached) {
            if (allowGoalNotReached) {
                return new Pair<Integer, Boolean>(CLIPPING_GOAL_NOT_REACHED, false);
            } else {
                throw new ReviewedStingException("Somehow the requested coordinate is not covered by the read. Too many deletions?");
            }
        }

        return new Pair<Integer, Boolean>(readBases, fallsInsideDeletion);
    }

    /**
     * Compares two SAMRecords only the basis on alignment start.  Note that
     * comparisons are performed ONLY on the basis of alignment start; any
     * two SAM records with the same alignment start will be considered equal.
     *
     * Unmapped alignments will all be considered equal.
     */

    @Requires({"read1 != null", "read2 != null"})
    public static int compareSAMRecords(GATKSAMRecord read1, GATKSAMRecord read2) {
        AlignmentStartComparator comp = new AlignmentStartComparator();
        return comp.compare(read1, read2);
    }

    /**
     * Is a base inside a read?
     *
     * @param read                the read to evaluate
     * @param referenceCoordinate the reference coordinate of the base to test
     * @return true if it is inside the read, false otherwise.
     */
    public static boolean isInsideRead(final GATKSAMRecord read, final int referenceCoordinate) {
        return referenceCoordinate >= read.getAlignmentStart() && referenceCoordinate <= read.getAlignmentEnd();
    }

    /**
     * Is this read all insertion?
     *
     * @param read
     * @return whether or not the only element in the cigar string is an Insertion
     */
    public static boolean readIsEntirelyInsertion(GATKSAMRecord read) {
        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() != CigarOperator.INSERTION)
                return false;
        }
        return true;
    }

    /**
     * Checks if a read starts with an insertion. It looks beyond Hard and Soft clips
     * if there are any.
     *
     * @param read
     * @return A pair with the answer (true/false) and the element or null if it doesn't exist
     */
    public static Pair<Boolean, CigarElement> readStartsWithInsertion(GATKSAMRecord read) {
        return readStartsWithInsertion(read.getCigar());
    }

    public static Pair<Boolean, CigarElement> readStartsWithInsertion(final Cigar cigar) {
        for (CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator() == CigarOperator.INSERTION)
                return new Pair<Boolean, CigarElement>(true, cigarElement);

            else if (cigarElement.getOperator() != CigarOperator.HARD_CLIP && cigarElement.getOperator() != CigarOperator.SOFT_CLIP)
                break;
        }
        return new Pair<Boolean, CigarElement>(false, null);
    }

    /**
     * Returns the coverage distribution of a list of reads within the desired region.
     *
     * See getCoverageDistributionOfRead for information on how the coverage is calculated.
     *
     * @param list          the list of reads covering the region
     * @param startLocation the first reference coordinate of the region (inclusive)
     * @param stopLocation  the last reference coordinate of the region (inclusive)
     * @return an array with the coverage of each position from startLocation to stopLocation
     */
    public static int [] getCoverageDistributionOfReads(List<GATKSAMRecord> list, int startLocation, int stopLocation) {
        int [] totalCoverage = new int[stopLocation - startLocation + 1];

        for (GATKSAMRecord read : list) {
            int [] readCoverage = getCoverageDistributionOfRead(read, startLocation, stopLocation);
            totalCoverage = MathUtils.addArrays(totalCoverage, readCoverage);
        }

        return totalCoverage;
    }

    /**
     * Returns the coverage distribution of a single read within the desired region.
     *
     * Note: This function counts DELETIONS as coverage (since the main purpose is to downsample
     * reads for variant regions, and deletions count as variants)
     *
     * @param read          the read to get the coverage distribution of
     * @param startLocation the first reference coordinate of the region (inclusive)
     * @param stopLocation  the last reference coordinate of the region (inclusive)
     * @return an array with the coverage of each position from startLocation to stopLocation
     */
    public static int [] getCoverageDistributionOfRead(GATKSAMRecord read, int startLocation, int stopLocation) {
        int [] coverage = new int[stopLocation - startLocation + 1];
        int refLocation = read.getSoftStart();
        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            switch (cigarElement.getOperator()) {
                case S:
                case M:
                case EQ:
                case N:
                case X:
                case D:
                    for (int i = 0; i < cigarElement.getLength(); i++) {
                        if (refLocation >= startLocation && refLocation <= stopLocation) {
                            int baseCount = read.isReducedRead() ? read.getReducedCount(refLocation - read.getSoftStart()) : 1;
                            coverage[refLocation - startLocation] += baseCount;   // this may be a reduced read, so add the proper number of bases
                        }
                        refLocation++;
                    }
                    break;

                case P:
                case I:
                case H:
                    break;
            }

            if (refLocation > stopLocation)
                break;
        }
        return coverage;
    }

    /**
     * Makes association maps for the reads and loci coverage as described below :
     *
     *  - First: locusToReadMap -- a HashMap that describes for each locus, which reads contribute to its coverage.
     *    Note: Locus is in reference coordinates.
     *    Example: Locus => {read1, read2, ..., readN}
     *
     *  - Second: readToLocusMap -- a HashMap that describes for each read what loci it contributes to the coverage.
     *    Note: Locus is a boolean array, indexed from 0 (= startLocation) to N (= stopLocation), with value==true meaning it contributes to the coverage.
     *    Example: Read => {true, true, false, ... false}
     *
     * @param readList      the list of reads to generate the association mappings
     * @param startLocation the first reference coordinate of the region (inclusive)
     * @param stopLocation  the last reference coordinate of the region (inclusive)
     * @return the two hashmaps described above
     */
    public static Pair<HashMap<Integer, HashSet<GATKSAMRecord>> , HashMap<GATKSAMRecord, Boolean[]>> getBothReadToLociMappings (List<GATKSAMRecord> readList, int startLocation, int stopLocation) {
        int arraySize = stopLocation - startLocation + 1;

        HashMap<Integer, HashSet<GATKSAMRecord>> locusToReadMap = new HashMap<Integer, HashSet<GATKSAMRecord>>(2*(stopLocation - startLocation + 1), 0.5f);
        HashMap<GATKSAMRecord, Boolean[]> readToLocusMap = new HashMap<GATKSAMRecord, Boolean[]>(2*readList.size(), 0.5f);

        for (int i = startLocation; i <= stopLocation; i++)
            locusToReadMap.put(i, new HashSet<GATKSAMRecord>()); // Initialize the locusToRead map with empty lists

        for (GATKSAMRecord read : readList) {
            readToLocusMap.put(read, new Boolean[arraySize]);       // Initialize the readToLocus map with empty arrays

            int [] readCoverage = getCoverageDistributionOfRead(read, startLocation, stopLocation);

            for (int i = 0; i < readCoverage.length; i++) {
                int refLocation = i + startLocation;
                if (readCoverage[i] > 0) {
                    // Update the hash for this locus
                    HashSet<GATKSAMRecord> readSet = locusToReadMap.get(refLocation);
                    readSet.add(read);

                    // Add this locus to the read hash
                    readToLocusMap.get(read)[refLocation - startLocation] = true;
                }
                else
                    // Update the boolean array with a 'no coverage' from this read to this locus
                    readToLocusMap.get(read)[refLocation-startLocation] = false;
            }
        }
        return new Pair<HashMap<Integer, HashSet<GATKSAMRecord>>, HashMap<GATKSAMRecord, Boolean[]>>(locusToReadMap, readToLocusMap);
    }

    /**
     * Create random read qualities
     *
     * @param length the length of the read
     * @return an array with randomized base qualities between 0 and 50
     */
    public static byte[] createRandomReadQuals(int length) {
        Random random = GenomeAnalysisEngine.getRandomGenerator();
        byte[] quals = new byte[length];
        for (int i = 0; i < length; i++)
            quals[i] = (byte) random.nextInt(50);
        return quals;
    }

    /**
     * Create random read qualities
     *
     * @param length  the length of the read
     * @param allowNs whether or not to allow N's in the read
     * @return an array with randomized bases (A-N) with equal probability
     */
    public static byte[] createRandomReadBases(int length, boolean allowNs) {
        Random random = GenomeAnalysisEngine.getRandomGenerator();
        int numberOfBases = allowNs ? 5 : 4;
        byte[] bases = new byte[length];
        for (int i = 0; i < length; i++) {
            switch (random.nextInt(numberOfBases)) {
                case 0:
                    bases[i] = 'A';
                    break;
                case 1:
                    bases[i] = 'C';
                    break;
                case 2:
                    bases[i] = 'G';
                    break;
                case 3:
                    bases[i] = 'T';
                    break;
                case 4:
                    bases[i] = 'N';
                    break;
                default:
                    throw new ReviewedStingException("Something went wrong, this is just impossible");
            }
        }
        return bases;
    }

    public static GATKSAMRecord createRandomRead(int length) {
        return createRandomRead(length, true);
    }

    public static GATKSAMRecord createRandomRead(int length, boolean allowNs) {
        byte[] quals = ReadUtils.createRandomReadQuals(length);
        byte[] bbases = ReadUtils.createRandomReadBases(length, allowNs);
        return ArtificialSAMUtils.createArtificialRead(bbases, quals, bbases.length + "M");
    }


    public static String prettyPrintSequenceRecords ( SAMSequenceDictionary sequenceDictionary ) {
        String[] sequenceRecordNames = new String[sequenceDictionary.size()];
        int sequenceRecordIndex = 0;
        for (SAMSequenceRecord sequenceRecord : sequenceDictionary.getSequences())
            sequenceRecordNames[sequenceRecordIndex++] = sequenceRecord.getSequenceName();
        return Arrays.deepToString(sequenceRecordNames);
    }

    /**
     * Calculates the reference coordinate for a read coordinate
     *
     * @param read   the read
     * @param offset the base in the read (coordinate in the read)
     * @return the reference coordinate correspondent to this base
     */
    public static long getReferenceCoordinateForReadCoordinate(GATKSAMRecord read, int offset) {
        if (offset > read.getReadLength()) 
            throw new ReviewedStingException(String.format(OFFSET_OUT_OF_BOUNDS_EXCEPTION, offset, read.getReadLength()));

        long location = read.getAlignmentStart();
        Iterator<CigarElement> cigarElementIterator = read.getCigar().getCigarElements().iterator();
        while (offset > 0 && cigarElementIterator.hasNext()) {
            CigarElement cigarElement = cigarElementIterator.next();
            long move = 0;
            if (cigarElement.getOperator().consumesReferenceBases())  
                move = (long) Math.min(cigarElement.getLength(), offset);
            location += move;
            offset -= move;
        }
        if (offset > 0 && !cigarElementIterator.hasNext()) 
            throw new ReviewedStingException(OFFSET_NOT_ZERO_EXCEPTION);

        return location;
    }

    /**
     * Creates a map with each event in the read (cigar operator) and the read coordinate where it happened.
     *
     * Example:
     *  D -> 2, 34, 75
     *  I -> 55
     *  S -> 0, 101
     *  H -> 101
     *
     * @param read the read
     * @return a map with the properties described above. See example
     */
    public static Map<CigarOperator, ArrayList<Integer>> getCigarOperatorForAllBases (GATKSAMRecord read) {
        Map<CigarOperator, ArrayList<Integer>> events = new HashMap<CigarOperator, ArrayList<Integer>>();

        int position = 0;
        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            CigarOperator op = cigarElement.getOperator();
            if (op.consumesReadBases()) {
                ArrayList<Integer> list = events.get(op);
                if (list == null) {
                    list = new ArrayList<Integer>();
                    events.put(op, list);
                }
                for (int i = position; i < cigarElement.getLength(); i++)
                    list.add(position++);
            }
            else {
                ArrayList<Integer> list = events.get(op);
                if (list == null) {
                    list = new ArrayList<Integer>();
                    events.put(op, list);
                }
                list.add(position);
            }
        }
        return events;
    }

}
