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
    private ReadUtils() { }

    private static int DEFAULT_ADAPTOR_SIZE = 100;

    public enum ClippingTail {
        LEFT_TAIL,
        RIGHT_TAIL
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
     * @return the reference coordinate for the adaptor boundary (effectively the first base IN the adaptor, closest to the read. NULL if the read is unmapped or the insert size cannot be determined (and is necessary for the calculation).
     */
    protected static Integer getAdaptorBoundary(final SAMRecord read) {
        if ( read.getReadUnmappedFlag() )
            return null;                                             // don't worry about unmapped pairs

        final int isize = Math.abs(read.getInferredInsertSize());    // the inferred insert size can be negative if the mate is mapped before the read (so we take the absolute value)
        int adaptorBoundary;                                         // the reference coordinate for the adaptor boundary (effectively the first base IN the adaptor, closest to the read)

        if ( read.getReadNegativeStrandFlag() )
            adaptorBoundary = read.getMateAlignmentStart() - 1;      // case 1 (see header)
        else if (isize > 0)
            adaptorBoundary = read.getAlignmentStart() + isize + 1;  // case 2 (see header)
        else
            return null;                                             // this is a case 2 where for some reason we cannot estimate the insert size

        return adaptorBoundary;
    }

    // return true if the read needs to be completely clipped
    private static GATKSAMRecord hardClipStartOfRead(GATKSAMRecord oldRec, int stopPosition) {

        if ( stopPosition >= oldRec.getAlignmentEnd() ) {
            // BAM representation issue -- we can't clip away all bases in a read, just leave it alone and let the filter deal with it
            //System.out.printf("Entire read needs to be clipped: %50s %n", read.getReadName());
            return null;
        }

        GATKSAMRecord read;
        try {
            read = (GATKSAMRecord)oldRec.clone();
        } catch (Exception e) {
            return null;
        }

        //System.out.printf("Clipping start of read: %50s start=%d adaptorEnd=%d isize=%d %n",
        //        read.getReadName(), read.getAlignmentStart(), stopPosition, read.getInferredInsertSize());

        Cigar oldCigar = read.getCigar();
        LinkedList<CigarElement> newCigarElements = new LinkedList<CigarElement>();
        int currentPos = read.getAlignmentStart();
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
        final byte[] bases = read.getReadBases();
        final byte[] quals = read.getBaseQualities();
        int newLength = bases.length - basesToClip;
        byte[] newBases = new byte[newLength];
        byte[] newQuals = new byte[newLength];
        System.arraycopy(bases, basesToClip, newBases, 0, newLength);
        System.arraycopy(quals, basesToClip, newQuals, 0, newLength);
        read.setReadBases(newBases);
        read.setBaseQualities(newQuals);

        // now add a CIGAR element for the clipped bases
        newCigarElements.addFirst(new CigarElement(basesToClip + basesAlreadyClipped, CigarOperator.H));
        Cigar newCigar = new Cigar(newCigarElements);
        read.setCigar(newCigar);

        // adjust the start accordingly
        read.setAlignmentStart(stopPosition + 1);

        return read;
    }

    private static GATKSAMRecord hardClipEndOfRead(GATKSAMRecord oldRec, int startPosition) {

        if ( startPosition <= oldRec.getAlignmentStart() ) {
            // BAM representation issue -- we can't clip away all bases in a read, just leave it alone and let the filter deal with it
            //System.out.printf("Entire read needs to be clipped: %50s %n", read.getReadName());
            return null;
        }

        GATKSAMRecord read;
        try {
            read = (GATKSAMRecord)oldRec.clone();
        } catch (Exception e) {
            return null;
        }

        //System.out.printf("Clipping end of read: %50s adaptorStart=%d end=%d isize=%d %n",
        //        read.getReadName(), startPosition, read.getAlignmentEnd(), read.getInferredInsertSize());

        Cigar oldCigar = read.getCigar();
        LinkedList<CigarElement> newCigarElements = new LinkedList<CigarElement>();
        int currentPos = read.getAlignmentStart();
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
        final byte[] bases = read.getReadBases();
        final byte[] quals = read.getBaseQualities();
        byte[] newBases = new byte[basesToKeep];
        byte[] newQuals = new byte[basesToKeep];
        System.arraycopy(bases, 0, newBases, 0, basesToKeep);
        System.arraycopy(quals, 0, newQuals, 0, basesToKeep);
        read.setReadBases(newBases);
        read.setBaseQualities(newQuals);

        // now add a CIGAR element for the clipped bases
        newCigarElements.add(new CigarElement((bases.length - basesToKeep) + basesAlreadyClipped, CigarOperator.H));
        Cigar newCigar = new Cigar(newCigarElements);
        read.setCigar(newCigar);

        // adjust the stop accordingly
        // read.setAlignmentEnd(startPosition - 1);

        return read;
    }

    /**
     * Hard clips away (i.e.g, removes from the read) bases that were previously soft clipped.
     *
     * @param read
     * @return
     */
    @Requires("read != null")
    @Ensures("result != null")
    public static GATKSAMRecord hardClipSoftClippedBases(GATKSAMRecord read) {
        List<CigarElement> cigarElts = read.getCigar().getCigarElements();

        if ( cigarElts.size() == 1 ) // can't be soft clipped, just return
            return read;

        int keepStart = 0, keepEnd = read.getReadLength() - 1;
        List<CigarElement> newCigarElements = new LinkedList<CigarElement>();

        for ( int i = 0; i < cigarElts.size(); i++ ) {
            CigarElement ce = cigarElts.get(i);
            int l = ce.getLength();
            switch ( ce.getOperator() ) {
                case S:
                    if ( i == 0 )
                        keepStart = l;
                    else
                        keepEnd = read.getReadLength() - l - 1;
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

        return hardClipBases(read, keepStart, keepEnd, mergedCigarElements);
    }

    /**
     * Hard clips out the bases in read, keeping the bases from keepStart to keepEnd, inclusive.  Note these
     * are offsets, so they are 0 based
     *
     * @param read
     * @param keepStart
     * @param keepEnd
     * @param newCigarElements
     * @return
     */
    @Requires({
            "read != null",
            "keepStart >= 0",
            "keepEnd < read.getReadLength()",
            "read.getReadUnmappedFlag() || newCigarElements != null"})
    @Ensures("result != null")
    public static GATKSAMRecord hardClipBases(GATKSAMRecord read, int keepStart, int keepEnd, List<CigarElement> newCigarElements) {
        int newLength = keepEnd - keepStart + 1;
        if ( newLength != read.getReadLength() ) {
            try {
                read = (GATKSAMRecord)read.clone();
                // copy over the unclipped bases
                final byte[] bases = read.getReadBases();
                final byte[] quals = read.getBaseQualities();
                byte[] newBases = new byte[newLength];
                byte[] newQuals = new byte[newLength];
                System.arraycopy(bases, keepStart, newBases, 0, newLength);
                System.arraycopy(quals, keepStart, newQuals, 0, newLength);
                read.setReadBases(newBases);
                read.setBaseQualities(newQuals);

                // now add a CIGAR element for the clipped bases, if the read isn't unmapped
                if ( ! read.getReadUnmappedFlag() ) {
                    Cigar newCigar = new Cigar(newCigarElements);
                    read.setCigar(newCigar);
                }
            } catch ( CloneNotSupportedException e ) {
                throw new ReviewedStingException("WTF, where did clone go?", e);
            }
        }

        return read;
    }

    public static GATKSAMRecord replaceSoftClipsWithMatches(GATKSAMRecord read) {
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



    /**
     *
     * @param read  original SAM record
     * @return a new read with adaptor sequence hard-clipped out or null if read is fully clipped
     */
    public static GATKSAMRecord hardClipAdaptorSequence(final GATKSAMRecord read) {
        final Integer adaptorBoundary = getAdaptorBoundary(read);

        if (adaptorBoundary == null || !isInsideRead(read, adaptorBoundary))
            return read;

        return (read.getReadNegativeStrandFlag()) ? hardClipStartOfRead(read, adaptorBoundary) : hardClipEndOfRead(read, adaptorBoundary);
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

    public static String readFlagsAsString(GATKSAMRecord read) {
        String flags = "";
        for (int flag : readFlagNames.keySet()) {
            if ((read.getFlags() & flag) != 0) {
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
    public final static List<GATKSAMRecord> coordinateSortReads(List<GATKSAMRecord> reads) {
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
    public static ReadAndIntervalOverlap getReadAndIntervalOverlapType(GATKSAMRecord read, GenomeLoc interval) {

        int sStart = getRefCoordSoftUnclippedStart(read);
        int sStop = getRefCoordSoftUnclippedEnd(read);
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

    @Ensures({"result >= read.getUnclippedStart()", "result <= read.getUnclippedEnd() || readIsEntirelyInsertion(read)"})
    public static int getRefCoordSoftUnclippedStart(GATKSAMRecord read) {
        int start = read.getUnclippedStart();
        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                start += cigarElement.getLength();
            else
                break;
        }
        return start;
    }

    @Ensures({"result >= read.getUnclippedStart()", "result <= read.getUnclippedEnd() || readIsEntirelyInsertion(read)"})
    public static int getRefCoordSoftUnclippedEnd(GATKSAMRecord read) {
        int stop = read.getUnclippedStart();

        if (readIsEntirelyInsertion(read))
            return stop;

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


    /**
     * Pre-processes the results of getReadCoordinateForReferenceCoordinate(GATKSAMRecord, int) in case it falls in
     * a deletion following the typical clipping needs. If clipping the left tail (beginning of the read) returns
     * the base prior to the deletion. If clipping the right tail (end of the read) returns the base after the
     * deletion.
     *
     * @param read
     * @param refCoord
     * @param tail
     * @return the read coordinate corresponding to the requested reference coordinate for clipping.
     */
    @Requires({"refCoord >= read.getUnclippedStart()", "refCoord <= read.getUnclippedEnd()"})
    @Ensures({"result >= 0", "result < read.getReadLength()"})
    public static int getReadCoordinateForReferenceCoordinate(GATKSAMRecord read, int refCoord, ClippingTail tail) {
        Pair<Integer, Boolean> result = getReadCoordinateForReferenceCoordinate(read, refCoord);
        int readCoord = result.getFirst();

        if (result.getSecond() && tail == ClippingTail.RIGHT_TAIL)
            readCoord++;

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
    @Requires({"refCoord >= getRefCoordSoftUnclippedStart(read)", "refCoord <= getRefCoordSoftUnclippedEnd(read)"})
    @Ensures({"result.getFirst() >= 0", "result.getFirst() < read.getReadLength()"})
    public static Pair<Integer, Boolean> getReadCoordinateForReferenceCoordinate(GATKSAMRecord read, int refCoord) {
        int readBases = 0;
        int refBases = 0;
        boolean fallsInsideDeletion = false;

        int goal = refCoord - getRefCoordSoftUnclippedStart(read);  // The goal is to move this many reference bases
        boolean goalReached = refBases == goal;

        Iterator<CigarElement> cigarElementIterator = read.getCigar().getCigarElements().iterator();
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
                boolean endsWithinCigar = shift <  cigarElement.getLength();

                // If it isn't, we need to check the next one. There should *ALWAYS* be a next one
                // since we checked if the goal coordinate is within the read length, so this is just a sanity check.
                if (!endsWithinCigar && !cigarElementIterator.hasNext())
                    throw new ReviewedStingException("Reference coordinate corresponds to a non-existent base in the read. This should never happen -- call Mauricio");

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
                        if (!cigarElementIterator.hasNext())
                            throw new ReviewedStingException("Reference coordinate corresponds to a non-existent base in the read. This should never happen -- call Mauricio");

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

            if (!goalReached)
                throw new ReviewedStingException("Somehow the requested coordinate is not covered by the read. Too many deletions?");
        

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
    @Ensures("result == 0 || result == 1 || result == -1")
    public static int compareSAMRecords(GATKSAMRecord read1, GATKSAMRecord read2) {
        AlignmentStartComparator comp = new AlignmentStartComparator();
        return comp.compare(read1, read2);
    }

    /**
     * Is a base inside a read?
     *
     * @param read the read to evaluate
     * @param referenceCoordinate the reference coordinate of the base to test
     * @return true if it is inside the read, false otherwise.
     */
    public static boolean isInsideRead(final GATKSAMRecord read, final int referenceCoordinate) {
        return referenceCoordinate >= read.getAlignmentStart() && referenceCoordinate <= read.getAlignmentEnd();
    }

    private static boolean readIsEntirelyInsertion(GATKSAMRecord read) {
        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() != CigarOperator.INSERTION)
                return false;
        }
        return true;
    }




}
