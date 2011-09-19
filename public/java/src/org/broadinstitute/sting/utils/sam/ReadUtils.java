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
import org.broadinstitute.sting.utils.GenomeLoc;
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
    public static final String REDUCED_READ_QUALITY_TAG = "RQ";

    private ReadUtils() { }

    public static SAMFileHeader copySAMFileHeader(SAMFileHeader toCopy) {
        SAMFileHeader copy = new SAMFileHeader();

        copy.setSortOrder(toCopy.getSortOrder());
        copy.setGroupOrder(toCopy.getGroupOrder());
        copy.setProgramRecords(toCopy.getProgramRecords());
        copy.setReadGroups(toCopy.getReadGroups());
        copy.setSequenceDictionary(toCopy.getSequenceDictionary());

        for (Map.Entry<String, String> e : toCopy.getAttributes())
            copy.setAttribute(e.getKey(), e.getValue());

        return copy;
    }

    public static SAMFileWriter createSAMFileWriterWithCompression(SAMFileHeader header, boolean presorted, String file, int compression) {
        if (file.endsWith(".bam"))
            return new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(file), compression);
        return new SAMFileWriterFactory().makeSAMOrBAMWriter(header, presorted, new File(file));
    }

    public static boolean isPlatformRead(SAMRecord read, String name) {
        SAMReadGroupRecord readGroup = read.getReadGroup();
        if (readGroup != null) {
            Object readPlatformAttr = readGroup.getAttribute("PL");
            if (readPlatformAttr != null)
                return readPlatformAttr.toString().toUpperCase().contains(name);
        }
        return false;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // utilities for detecting overlapping reads
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Detects read pairs where the reads are so long relative to the over fragment size that they are
     * reading into each other's adaptors.
     *
     * Normally, fragments are sufficiently far apart that reads aren't reading into each other.
     *
     * |-------------------->                                   first read
     *                                 <--------------------|   second read
     *
     * Sometimes, mostly due to lab errors or constraints, fragment library are made too short relative to the
     * length of the reads.  For example, it's possible to have 76bp PE reads with 125 bp inserts, so that ~25 bp of each
     * read overlaps with its mate.
     *
     * |--------OOOOOOOOOOOO>               first read
     *         <OOOOOOOOOOOO------------|   second read
     *
     * This filter deals with the situation where the fragment is so small that the each read actually reads into the
     * adaptor sequence of its mate, generating mismatches at both ends of the read:
     *
     *              |----------------XXXX>      first read
     *         <XXXX----------------|           second read
     *
     * The code below returns NOT_OVERLAPPING for the first case, IN_ADAPTOR for the last case, and OVERLAPPING
     * given a read and a reference aligned base position.
     *
     * @author depristo
     * @version 0.1
     */

    public enum OverlapType { NOT_OVERLAPPING, IN_ADAPTOR}

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
    public enum ReadAndIntervalOverlap {NO_OVERLAP_CONTIG, NO_OVERLAP_LEFT, NO_OVERLAP_RIGHT, OVERLAP_LEFT, OVERLAP_RIGHT, OVERLAP_LEFT_AND_RIGHT, OVERLAP_CONTAINED}

    /**
     * God, there's a huge information asymmetry in SAM format:
     *
     *      s1                      e1
     *      |-----------------------> [record in hand]
     *  s2
     *  <-----------------------|
     *
     * s1, e1, and s2 are all in the record.  From isize we can can compute e2 as s1 + isize + 1
     *
     *      s2
     *      |----------------------->
     *  s1                      e1
     *  <-----------------------|     [record in hand]
     *
     * Here we cannot calculate e2 since the record carries s2 and e1 + isize is s2 now!
     *
     * This makes the following code a little nasty, since we can only detect if a base is in the adaptor, but not
     * if it overlaps the read.
     *
     * @param rec
     * @param basePos
     * @param adaptorLength
     * @return
     */
    public static OverlapType readPairBaseOverlapType(final SAMRecord rec, long basePos, final int adaptorLength) {
        OverlapType state = OverlapType.NOT_OVERLAPPING;

        Pair<Integer, Integer> adaptorBoundaries = getAdaptorBoundaries(rec, adaptorLength);

        if ( adaptorBoundaries != null ) { // we're not an unmapped pair -- cannot filter out

            boolean inAdapator = basePos >= adaptorBoundaries.first && basePos <= adaptorBoundaries.second;

            if ( inAdapator ) { 
                state = OverlapType.IN_ADAPTOR;
                //System.out.printf("baseOverlapState: %50s negStrand=%b base=%d start=%d stop=%d, adaptorStart=%d adaptorEnd=%d isize=%d => %s%n",
                //        rec.getReadName(), rec.getReadNegativeStrandFlag(), basePos, rec.getAlignmentStart(), rec.getAlignmentEnd(), adaptorBoundaries.first, adaptorBoundaries.second, rec.getInferredInsertSize(), state);
            }
        }

        return state;
    }

    private static Pair<Integer, Integer> getAdaptorBoundaries(SAMRecord rec, int adaptorLength) {
        int isize = rec.getInferredInsertSize();
        if ( isize == 0 )
            return null; // don't worry about unmapped pairs

        int adaptorStart, adaptorEnd;

        if ( rec.getReadNegativeStrandFlag() ) {
            // we are on the negative strand, so our mate is on the positive strand
            int mateStart = rec.getMateAlignmentStart();
            adaptorStart = mateStart - adaptorLength - 1;
            adaptorEnd = mateStart - 1;
        } else {
            // we are on the positive strand, so our mate is on the negative strand
            int mateEnd = rec.getAlignmentStart() + isize - 1;
            adaptorStart = mateEnd + 1;
            adaptorEnd = mateEnd + adaptorLength;
        }

        return new Pair<Integer, Integer>(adaptorStart, adaptorEnd);
    }

    /**
     *
     * @param rec  original SAM record
     * @param adaptorLength  length of adaptor sequence
     * @return a new read with adaptor sequence hard-clipped out or null if read is fully clipped
     */
    public static GATKSAMRecord hardClipAdaptorSequence(final SAMRecord rec, int adaptorLength) {

        Pair<Integer, Integer> adaptorBoundaries = getAdaptorBoundaries(rec, adaptorLength);
        GATKSAMRecord result = (GATKSAMRecord)rec;

        if ( adaptorBoundaries != null ) {
            if ( rec.getReadNegativeStrandFlag() && adaptorBoundaries.second >= rec.getAlignmentStart() && adaptorBoundaries.first < rec.getAlignmentEnd() )
                result = hardClipStartOfRead(rec, adaptorBoundaries.second);
            else if ( !rec.getReadNegativeStrandFlag() && adaptorBoundaries.first <= rec.getAlignmentEnd() )
                result = hardClipEndOfRead(rec, adaptorBoundaries.first);
        }

        return result;
    }

    // return true if the read needs to be completely clipped
    private static GATKSAMRecord hardClipStartOfRead(SAMRecord oldRec, int stopPosition) {

        if ( stopPosition >= oldRec.getAlignmentEnd() ) {
            // BAM representation issue -- we can't clip away all bases in a read, just leave it alone and let the filter deal with it
            //System.out.printf("Entire read needs to be clipped: %50s %n", rec.getReadName());
            return null;
        }

        GATKSAMRecord rec;
        try {
            rec = (GATKSAMRecord)oldRec.clone();
        } catch (Exception e) {
            return null;
        }

        //System.out.printf("Clipping start of read: %50s start=%d adaptorEnd=%d isize=%d %n",
        //        rec.getReadName(), rec.getAlignmentStart(), stopPosition, rec.getInferredInsertSize());

        Cigar oldCigar = rec.getCigar();
        LinkedList<CigarElement> newCigarElements = new LinkedList<CigarElement>();
        int currentPos = rec.getAlignmentStart();
        int basesToClip = 0;
        int basesAlreadyClipped = 0;

        for ( CigarElement ce : oldCigar.getCigarElements() ) {

            if ( currentPos > stopPosition) {
                newCigarElements.add(ce);
                continue;
            }

            int elementLength = ce.getLength();
            switch ( ce.getOperator() ) {
                case M:
                    for (int i = 0; i < elementLength; i++, currentPos++, basesToClip++) {
                        if ( currentPos > stopPosition ) {
                            newCigarElements.add(new CigarElement(elementLength - i, CigarOperator.M));
                            break;
                        }
                    }
                    break;
                case I:
                case S:
                    basesToClip += elementLength;
                    break;
                case D:
                case N:
                    currentPos += elementLength;
                    break;
                case H:
                    basesAlreadyClipped += elementLength;
                case P:
                    break;
                default: throw new ReviewedStingException("The " + ce.getOperator() + " cigar element is not currently supported");
            }

        }

        // copy over the unclipped bases
        final byte[] bases = rec.getReadBases();
        final byte[] quals = rec.getBaseQualities();
        int newLength = bases.length - basesToClip;
        byte[] newBases = new byte[newLength];
        byte[] newQuals = new byte[newLength];
        System.arraycopy(bases, basesToClip, newBases, 0, newLength);
        System.arraycopy(quals, basesToClip, newQuals, 0, newLength);
        rec.setReadBases(newBases);
        rec.setBaseQualities(newQuals);

        // now add a CIGAR element for the clipped bases
        newCigarElements.addFirst(new CigarElement(basesToClip + basesAlreadyClipped, CigarOperator.H));
        Cigar newCigar = new Cigar(newCigarElements);
        rec.setCigar(newCigar);

        // adjust the start accordingly
        rec.setAlignmentStart(stopPosition + 1);

        return rec;
    }

    private static GATKSAMRecord hardClipEndOfRead(SAMRecord oldRec, int startPosition) {

        if ( startPosition <= oldRec.getAlignmentStart() ) {
            // BAM representation issue -- we can't clip away all bases in a read, just leave it alone and let the filter deal with it
            //System.out.printf("Entire read needs to be clipped: %50s %n", rec.getReadName());
            return null;
        }

        GATKSAMRecord rec;
        try {
            rec = (GATKSAMRecord)oldRec.clone();
        } catch (Exception e) {
            return null;
        }

        //System.out.printf("Clipping end of read: %50s adaptorStart=%d end=%d isize=%d %n",
        //        rec.getReadName(), startPosition, rec.getAlignmentEnd(), rec.getInferredInsertSize());

        Cigar oldCigar = rec.getCigar();
        LinkedList<CigarElement> newCigarElements = new LinkedList<CigarElement>();
        int currentPos = rec.getAlignmentStart();
        int basesToKeep = 0;
        int basesAlreadyClipped = 0;

        for ( CigarElement ce : oldCigar.getCigarElements() ) {

            int elementLength = ce.getLength();

            if ( currentPos >= startPosition ) {
                if ( ce.getOperator() == CigarOperator.H )
                    basesAlreadyClipped += elementLength;
                continue;
            }

            switch ( ce.getOperator() ) {
                case M:
                    for (int i = 0; i < elementLength; i++, currentPos++, basesToKeep++) {
                        if ( currentPos == startPosition ) {
                            newCigarElements.add(new CigarElement(i, CigarOperator.M));
                            break;
                        }
                    }

                    if ( currentPos != startPosition )
                        newCigarElements.add(ce);
                    break;
                case I:
                case S:
                    newCigarElements.add(ce);
                    basesToKeep += elementLength;
                    break;
                case D:
                case N:
                    newCigarElements.add(ce);
                    currentPos += elementLength;
                    break;
                case H:
                case P:
                    newCigarElements.add(ce);
                    break;
                default: throw new ReviewedStingException("The " + ce.getOperator() + " cigar element is not currently supported");
            }

        }

        // copy over the unclipped bases
        final byte[] bases = rec.getReadBases();
        final byte[] quals = rec.getBaseQualities();
        byte[] newBases = new byte[basesToKeep];
        byte[] newQuals = new byte[basesToKeep];
        System.arraycopy(bases, 0, newBases, 0, basesToKeep);
        System.arraycopy(quals, 0, newQuals, 0, basesToKeep);
        rec.setReadBases(newBases);
        rec.setBaseQualities(newQuals);

        // now add a CIGAR element for the clipped bases
        newCigarElements.add(new CigarElement((bases.length - basesToKeep) + basesAlreadyClipped, CigarOperator.H));
        Cigar newCigar = new Cigar(newCigarElements);
        rec.setCigar(newCigar);

        // adjust the stop accordingly
        // rec.setAlignmentEnd(startPosition - 1);

        return rec;
    }

    /**
     * Hard clips away (i.e.g, removes from the read) bases that were previously soft clipped.
     *
     * @param rec
     * @return
     */
    @Requires("rec != null")
    @Ensures("result != null")
    public static SAMRecord hardClipSoftClippedBases(SAMRecord rec) {
        List<CigarElement> cigarElts = rec.getCigar().getCigarElements();

        if ( cigarElts.size() == 1 ) // can't be soft clipped, just return
            return rec;

        int keepStart = 0, keepEnd = rec.getReadLength() - 1;
        List<CigarElement> newCigarElements = new LinkedList<CigarElement>();

        for ( int i = 0; i < cigarElts.size(); i++ ) {
            CigarElement ce = cigarElts.get(i);
            int l = ce.getLength();
            switch ( ce.getOperator() ) {
                case S:
                    if ( i == 0 )
                        keepStart = l;
                    else
                        keepEnd = rec.getReadLength() - l - 1;
                    newCigarElements.add(new CigarElement(l, CigarOperator.HARD_CLIP));
                    break;

                default:
                    newCigarElements.add(ce);
                    break;
            }
        }

        // Merges tandem cigar elements like 5H10H or 2S5S to 15H or 7S
        // this will happen if you soft clip a read that has been hard clipped before
        // like: 5H20S => 5H20H
        List<CigarElement> mergedCigarElements = new LinkedList<CigarElement>();
        Iterator<CigarElement> cigarElementIterator = newCigarElements.iterator();
        CigarOperator currentOperator = null;
        int currentOperatorLength = 0;
        while (cigarElementIterator.hasNext()) {
            CigarElement cigarElement = cigarElementIterator.next();
            if (currentOperator != cigarElement.getOperator()) {
                if (currentOperator != null)
                    mergedCigarElements.add(new CigarElement(currentOperatorLength, currentOperator));
                currentOperator = cigarElement.getOperator();
                currentOperatorLength = cigarElement.getLength();
            }
            else
                currentOperatorLength += cigarElement.getLength();
        }
        mergedCigarElements.add(new CigarElement(currentOperatorLength, currentOperator));

        return hardClipBases(rec, keepStart, keepEnd, mergedCigarElements);
    }

    /**
     * Hard clips out the bases in rec, keeping the bases from keepStart to keepEnd, inclusive.  Note these
     * are offsets, so they are 0 based
     *
     * @param rec
     * @param keepStart
     * @param keepEnd
     * @param newCigarElements
     * @return
     */
    @Requires({
            "rec != null",
            "keepStart >= 0",
            "keepEnd < rec.getReadLength()",
            "rec.getReadUnmappedFlag() || newCigarElements != null"})
    @Ensures("result != null")
    public static SAMRecord hardClipBases(SAMRecord rec, int keepStart, int keepEnd, List<CigarElement> newCigarElements) {
        int newLength = keepEnd - keepStart + 1;
        if ( newLength != rec.getReadLength() ) {
            try {
                rec = SimplifyingSAMFileWriter.simplifyRead((SAMRecord)rec.clone());
                // copy over the unclipped bases
                final byte[] bases = rec.getReadBases();
                final byte[] quals = rec.getBaseQualities();
                byte[] newBases = new byte[newLength];
                byte[] newQuals = new byte[newLength];
                System.arraycopy(bases, keepStart, newBases, 0, newLength);
                System.arraycopy(quals, keepStart, newQuals, 0, newLength);
                rec.setReadBases(newBases);
                rec.setBaseQualities(newQuals);

                // now add a CIGAR element for the clipped bases, if the read isn't unmapped
                if ( ! rec.getReadUnmappedFlag() ) {
                    Cigar newCigar = new Cigar(newCigarElements);
                    rec.setCigar(newCigar);
                }
            } catch ( CloneNotSupportedException e ) {
                throw new ReviewedStingException("WTF, where did clone go?", e);
            }
        }

        return rec;
    }

    public static SAMRecord replaceSoftClipsWithMatches(SAMRecord read) {
        List<CigarElement> newCigarElements = new ArrayList<CigarElement>();

        for ( CigarElement ce : read.getCigar().getCigarElements() ) {
            if ( ce.getOperator() == CigarOperator.SOFT_CLIP )
                newCigarElements.add(new CigarElement(ce.getLength(), CigarOperator.MATCH_OR_MISMATCH));
            else
                newCigarElements.add(ce);
        }

        if ( newCigarElements.size() > 1 ) { //
            CigarElement first = newCigarElements.get(0);
            CigarElement second = newCigarElements.get(1);
            if ( first.getOperator() == CigarOperator.MATCH_OR_MISMATCH && second.getOperator() == CigarOperator.MATCH_OR_MISMATCH ) {
                newCigarElements.set(0, new CigarElement(first.getLength() + second.getLength(), CigarOperator.MATCH_OR_MISMATCH));
                newCigarElements.remove(1);
            }
        }

        if ( newCigarElements.size() > 1 ) { //
            CigarElement penult = newCigarElements.get(newCigarElements.size()-2);
            CigarElement last = newCigarElements.get(newCigarElements.size()-1);
            if ( penult.getOperator() == CigarOperator.MATCH_OR_MISMATCH && penult.getOperator() == CigarOperator.MATCH_OR_MISMATCH ) {
                newCigarElements.set(newCigarElements.size()-2, new CigarElement(penult.getLength() + last.getLength(), CigarOperator.MATCH_OR_MISMATCH));
                newCigarElements.remove(newCigarElements.size()-1);
            }
        }

        read.setCigar(new Cigar(newCigarElements));
        return read;
    }


    private static int DEFAULT_ADAPTOR_SIZE = 100;

    /**
     *
     * @param rec  original SAM record
     * @return a new read with adaptor sequence hard-clipped out or null if read is fully clipped
     */
    public static GATKSAMRecord hardClipAdaptorSequence(final SAMRecord rec) {
        return hardClipAdaptorSequence(rec, DEFAULT_ADAPTOR_SIZE);
    }

    public static OverlapType readPairBaseOverlapType(final SAMRecord rec, long basePos) {
        return readPairBaseOverlapType(rec, basePos, DEFAULT_ADAPTOR_SIZE);
    }

    public static boolean is454Read(SAMRecord read) {
        return isPlatformRead(read, "454");
    }

    public static boolean isSOLiDRead(SAMRecord read) {
        return isPlatformRead(read, "SOLID");
    }

    public static boolean isSLXRead(SAMRecord read) {
        return isPlatformRead(read, "ILLUMINA");
    }

    private static final Map<Integer, String> readFlagNames
            = new HashMap<Integer, String>();

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

    public static String readFlagsAsString(SAMRecord rec) {
        String flags = "";
        for (int flag : readFlagNames.keySet()) {
            if ((rec.getFlags() & flag) != 0) {
                flags += readFlagNames.get(flag) + " ";
            }
        }
        return flags;
    }

    /**
     * Returns the collections of reads sorted in coordinate order, according to the order defined
     * in the reads themselves
     *
     * @param reads
     * @return
     */
    public final static List<SAMRecord> coordinateSortReads(List<SAMRecord> reads) {
        final SAMRecordComparator comparer = new SAMRecordCoordinateComparator();
        Collections.sort(reads, comparer);
        return reads;
    }

    public final static int getFirstInsertionOffset(SAMRecord read) {
        CigarElement e = read.getCigar().getCigarElement(0);
        if ( e.getOperator() == CigarOperator.I )
            return e.getLength();
        else
            return 0;
    }

    public final static int getLastInsertionOffset(SAMRecord read) {
        CigarElement e = read.getCigar().getCigarElement(read.getCigarLength()-1);
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
    public static ReadAndIntervalOverlap getReadAndIntervalOverlapType(SAMRecord read, GenomeLoc interval) {

        int start = getRefCoordSoftUnclippedStart(read);
        int stop = getRefCoordSoftUnclippedEnd(read);

        if ( !read.getReferenceName().equals(interval.getContig()) )
            return ReadAndIntervalOverlap.NO_OVERLAP_CONTIG;

        else if  ( stop < interval.getStart() )
            return ReadAndIntervalOverlap.NO_OVERLAP_LEFT;

        else if ( start > interval.getStop() )
            return ReadAndIntervalOverlap.NO_OVERLAP_RIGHT;

        else if ( (start >= interval.getStart()) &&
                  (stop <= interval.getStop()) )
            return ReadAndIntervalOverlap.OVERLAP_CONTAINED;

        else if ( (start < interval.getStart()) &&
                  (stop > interval.getStop()) )
            return ReadAndIntervalOverlap.OVERLAP_LEFT_AND_RIGHT;

        else if ( (start < interval.getStart()) )
            return ReadAndIntervalOverlap.OVERLAP_LEFT;

        else
            return ReadAndIntervalOverlap.OVERLAP_RIGHT;
    }

    @Ensures({"(result >= read.getUnclippedStart() && result <= read.getUnclippedEnd()) || readIsEntirelyInsertion(read)"})
    public static int getRefCoordSoftUnclippedStart(SAMRecord read) {
        int start = read.getUnclippedStart();
        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                start += cigarElement.getLength();
            else
                break;
        }
        return start;
    }

    @Ensures({"(result >= read.getUnclippedStart() && result <= read.getUnclippedEnd()) || readIsEntirelyInsertion(read)"})
    public static int getRefCoordSoftUnclippedEnd(SAMRecord read) {
        int stop = read.getUnclippedStart();
        int shift = 0;
        CigarOperator lastOperator = null;
        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            stop += shift;
            lastOperator = cigarElement.getOperator();
            if (cigarElement.getOperator().consumesReferenceBases() || cigarElement.getOperator() == CigarOperator.SOFT_CLIP || cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                shift = cigarElement.getLength();
            else
                shift = 0;
        }
        return (lastOperator == CigarOperator.HARD_CLIP) ? stop-1 : stop+shift-1 ;
    }

    private static boolean readIsEntirelyInsertion(SAMRecord read) {
        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() != CigarOperator.INSERTION)
                return false;
        }
        return true;
    }

    /**
     * Looks for a read coordinate that corresponds to the reference coordinate in the soft clipped region before
     * the alignment start of the read.
     *
     * @param read
     * @param refCoord
     * @return the corresponding read coordinate or -1 if it failed to find it (it has been hard clipped before)
     */
    @Requires({"refCoord >= read.getUnclippedStart()", "refCoord < read.getAlignmentStart()"})
    private static int getReadCoordinateForReferenceCoordinateBeforeAlignmentStart(SAMRecord read, int refCoord) {
        if (getRefCoordSoftUnclippedStart(read) <= refCoord)
            return refCoord - getRefCoordSoftUnclippedStart(read) + 1;
        return -1;
    }


    /**
     * Looks for a read coordinate that corresponds to the reference coordinate in the soft clipped region after
     * the alignment end of the read.
     *
     * @param read
     * @param refCoord
     * @return the corresponding read coordinate or -1 if it failed to find it (it has been hard clipped before)
     */
    @Requires({"refCoord <= read.getUnclippedEnd()", "refCoord > read.getAlignmentEnd()"})
    private static int getReadCoordinateForReferenceCoordinateBeforeAlignmentEnd(SAMRecord read, int refCoord) {
        if (getRefCoordSoftUnclippedEnd(read) >= refCoord)
            return refCoord - getRefCoordSoftUnclippedStart(read) + 1;
        return -1;
    }


    @Requires({"refCoord >= read.getUnclippedStart()", "refCoord <= read.getUnclippedEnd()"})
    @Ensures({"result >= 0", "result < read.getReadLength()"})
    public static int getReadCoordinateForReferenceCoordinate(SAMRecord read, int refCoord) {
        int readBases = 0;
        int refBases = 0;

        if (refCoord < read.getAlignmentStart()) {
            readBases = getReadCoordinateForReferenceCoordinateBeforeAlignmentStart(read, refCoord);
            if (readBases < 0)
                throw new ReviewedStingException("Requested a coordinate in a hard clipped area of the read. No equivalent read coordinate.");
        }
        else if (refCoord > read.getAlignmentEnd()) {
            readBases = getReadCoordinateForReferenceCoordinateBeforeAlignmentEnd(read, refCoord);
            if (readBases < 0)
                throw new ReviewedStingException("Requested a coordinate in a hard clipped area of the read. No equivalent read coordinate.");
        }
        else {
            int goal = refCoord - read.getAlignmentStart();  // The goal is to move this many reference bases
            boolean goalReached = refBases == goal;

            Iterator<CigarElement> cigarElementIterator = read.getCigar().getCigarElements().iterator();
            while (!goalReached && cigarElementIterator.hasNext()) {
                CigarElement cigarElement = cigarElementIterator.next();
                int shift = 0;

                if (cigarElement.getOperator().consumesReferenceBases()) {
                    if (refBases + cigarElement.getLength() < goal) {
                        shift = cigarElement.getLength();
                    }
                    else {
                        shift = goal - refBases;
                    }
                    refBases += shift;
                }
                goalReached = refBases == goal;

                if (cigarElement.getOperator().consumesReadBases()) {
                    readBases += goalReached ? shift : cigarElement.getLength();
                }
            }

            if (!goalReached)
                throw new ReviewedStingException("Somehow the requested coordinate is not covered by the read. Too many deletions?");
        }

        return readBases;
    }


}
