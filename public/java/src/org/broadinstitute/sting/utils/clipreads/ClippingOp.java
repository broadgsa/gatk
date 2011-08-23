package org.broadinstitute.sting.utils.clipreads;

import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.apache.poi.hssf.record.PageBreakRecord;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.PipedOutputStream;
import java.lang.reflect.Array;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

/**
 * Represents a clip on a read.  It has a type (see the enum) along with a start and stop in the bases
 * of the read, plus an option extraInfo (useful for carrying info where needed).
 * <p/>
 * Also holds the critical apply function that actually execute the clipping operation on a provided read,
 * according to the wishes of the supplid ClippingAlgorithm enum.
 */
public class ClippingOp {
    public final int start, stop; // inclusive

    public ClippingOp(int start, int stop) {
        this.start = start;
        this.stop = stop;
    }


    public int getLength() {
        return stop - start + 1;
    }

    /**
     * Clips the bases in read according to this operation's start and stop.  Uses the clipping
     * representation used is the one provided by algorithm argument.
     *
     * @param algorithm
     * @param read
     */
    public SAMRecord apply(ClippingRepresentation algorithm, SAMRecord read) {
        byte[] quals = read.getBaseQualities();
        byte[] bases = read.getReadBases();

        switch (algorithm) {
            // important note:
            //   it's not safe to call read.getReadBases()[i] = 'N' or read.getBaseQualities()[i] = 0
            //   because you're not guaranteed to get a pointer to the actual array of bytes in the SAMRecord
            case WRITE_NS:
                for (int i = start; i <= stop; i++)
                    bases[i] = 'N';
                read.setReadBases(bases);
                break;
            case WRITE_Q0S:
                for (int i = start; i <= stop; i++)
                    quals[i] = 0;
                read.setBaseQualities(quals);
                break;
            case WRITE_NS_Q0S:
                for (int i = start; i <= stop; i++) {
                    bases[i] = 'N';
                    quals[i] = 0;
                }
                read.setReadBases(bases);
                read.setBaseQualities(quals);
                break;
            case HARDCLIP_BASES:
                read = hardClip(read, start, stop);
                break;

            case SOFTCLIP_BASES:
                if ( read.getReadUnmappedFlag() ) {
                    // we can't process unmapped reads
                    throw new UserException("Read Clipper cannot soft clip unmapped reads");
                }

                //System.out.printf("%d %d %d%n", stop, start, read.getReadLength());
                int myStop = stop;
                if ( (stop + 1 - start) == read.getReadLength() ) {
                    // BAM representation issue -- we can't SOFTCLIP away all bases in a read, just leave it alone
                    //Walker.logger.info(String.format("Warning, read %s has all bases clip but this can't be represented with SOFTCLIP_BASES, just leaving it alone", read.getReadName()));
                    //break;
                    myStop--; // just decrement stop
                }

                if ( start > 0 && myStop != read.getReadLength() - 1 )
                    throw new RuntimeException(String.format("Cannot apply soft clipping operator to the middle of a read: %s to be clipped at %d-%d", read.getReadName(), start, myStop));

                Cigar oldCigar = read.getCigar();

                int scLeft = 0, scRight = read.getReadLength();
                if ( start == 0 )
                    scLeft = myStop + 1;
                else
                    scRight = start;

                Cigar newCigar = softClip(oldCigar, scLeft, scRight);
                read.setCigar(newCigar);

                int newClippedStart = getNewAlignmentStartOffset(newCigar, oldCigar);
                int newStart = read.getAlignmentStart() + newClippedStart;
                read.setAlignmentStart(newStart);

                break;

            default:
                throw new IllegalStateException("Unexpected Clipping operator type " + algorithm);
        }

        return read;
    }

    /**
     * Given a cigar string, get the number of bases hard or soft clipped at the start
     */
    private int getNewAlignmentStartOffset(final Cigar __cigar, final Cigar __oldCigar) {
        int num = 0;
        for (CigarElement e : __cigar.getCigarElements()) {
            if (!e.getOperator().consumesReferenceBases()) {
                if (e.getOperator().consumesReadBases()) {
                    num += e.getLength();
                }
            } else {
                break;
            }
        }

        int oldNum = 0;
        int curReadCounter = 0;

        for (CigarElement e : __oldCigar.getCigarElements()) {
            int curRefLength = e.getLength();
            int curReadLength = e.getLength();
            if (!e.getOperator().consumesReadBases()) {
                curReadLength = 0;
            }

            boolean truncated = false;
            if (curReadCounter + curReadLength > num) {
                curReadLength = num - curReadCounter;
                curRefLength = num - curReadCounter;
                truncated = true;
            }

            if (!e.getOperator().consumesReferenceBases()) {
                curRefLength = 0;
            }

            curReadCounter += curReadLength;
            oldNum += curRefLength;

            if (curReadCounter > num || truncated) {
                break;
            }
        }

        return oldNum;
    }

    /**
     * Given a cigar string, soft clip up to startClipEnd and soft clip starting at endClipBegin
     */
    private Cigar softClip(final Cigar __cigar, final int __startClipEnd, final int __endClipBegin) {
        if (__endClipBegin <= __startClipEnd) {
            //whole thing should be soft clipped
            int cigarLength = 0;
            for (CigarElement e : __cigar.getCigarElements()) {
                cigarLength += e.getLength();
            }

            Cigar newCigar = new Cigar();
            newCigar.add(new CigarElement(cigarLength, CigarOperator.SOFT_CLIP));
            assert newCigar.isValid(null, -1) == null;
            return newCigar;
        }

        int curLength = 0;
        Vector<CigarElement> newElements = new Vector<CigarElement>();
        for (CigarElement curElem : __cigar.getCigarElements()) {
            if (!curElem.getOperator().consumesReadBases()) {
                if (curElem.getOperator() == CigarOperator.HARD_CLIP || curLength > __startClipEnd && curLength < __endClipBegin) {
                    newElements.add(new CigarElement(curElem.getLength(), curElem.getOperator()));
                }
                continue;
            }

            int s = curLength;
            int e = curLength + curElem.getLength();
            if (e <= __startClipEnd || s >= __endClipBegin) {
                //must turn this entire thing into a clip
                newElements.add(new CigarElement(curElem.getLength(), CigarOperator.SOFT_CLIP));
            } else if (s >= __startClipEnd && e <= __endClipBegin) {
                //same thing
                newElements.add(new CigarElement(curElem.getLength(), curElem.getOperator()));
            } else {
                //we are clipping in the middle of this guy
                CigarElement newStart = null;
                CigarElement newMid = null;
                CigarElement newEnd = null;

                int midLength = curElem.getLength();
                if (s < __startClipEnd) {
                    newStart = new CigarElement(__startClipEnd - s, CigarOperator.SOFT_CLIP);
                    midLength -= newStart.getLength();
                }

                if (e > __endClipBegin) {
                    newEnd = new CigarElement(e - __endClipBegin, CigarOperator.SOFT_CLIP);
                    midLength -= newEnd.getLength();
                }
                assert midLength >= 0;
                if (midLength > 0) {
                    newMid = new CigarElement(midLength, curElem.getOperator());
                }
                if (newStart != null) {
                    newElements.add(newStart);
                }
                if (newMid != null) {
                    newElements.add(newMid);
                }
                if (newEnd != null) {
                    newElements.add(newEnd);
                }
            }
            curLength += curElem.getLength();
        }

        Vector<CigarElement> finalNewElements = new Vector<CigarElement>();
        CigarElement lastElement = null;
        for (CigarElement elem : newElements) {
            if (lastElement == null || lastElement.getOperator() != elem.getOperator()) {
                if (lastElement != null) {
                    finalNewElements.add(lastElement);
                }
                lastElement = elem;
            } else {
                lastElement = new CigarElement(lastElement.getLength() + elem.getLength(), lastElement.getOperator());
            }
        }
        if (lastElement != null) {
            finalNewElements.add(lastElement);
        }

        Cigar newCigar = new Cigar(finalNewElements);
        assert newCigar.isValid(null, -1) == null;
        return newCigar;
    }

    @Requires({"start <= stop", "start == 0 || stop == read.getReadLength() - 1", "!read.getReadUnmappedFlag()"})
    private SAMRecord hardClip (SAMRecord read, int start, int stop) {
        if (start == 0 && stop == read.getReadLength() -1)
            return new SAMRecord(read.getHeader());

        int newLength = read.getReadLength() - (stop - start + 1);
        byte [] newBases = new byte[newLength];
        byte [] newQuals = new byte[newLength];
        int copyStart = (start == 0) ? stop + 1 : 0;

        System.arraycopy(read.getReadBases(), copyStart, newBases, 0, newLength);
        System.arraycopy(read.getBaseQualities(), copyStart, newQuals, 0, newLength);        Cigar newCigar = hardClipCigar(read.getCigar(), start, stop);

        SAMRecord hardClippedRead;
        try {
            hardClippedRead = (SAMRecord) read.clone();
        } catch (CloneNotSupportedException e) {
            throw new ReviewedStingException("Where did the clone go?");
        }

        hardClippedRead.setBaseQualities(newQuals);
        hardClippedRead.setReadBases(newBases);
        hardClippedRead.setCigar(newCigar);
        if (start == 0)
            hardClippedRead.setAlignmentStart(read.getAlignmentStart() + calculateAlignmentStartShift(read.getCigar(), newCigar));

        return hardClippedRead;

    }

    @Requires({"!cigar.isEmpty()"})
    private Cigar hardClipCigar (Cigar cigar, int start, int stop) {
        Cigar newCigar = new Cigar();
        int index = 0;
        int totalHardClipCount = stop - start + 1;

        // hard clip the beginning of the cigar string
        if (start == 0) {
            Iterator<CigarElement> cigarElementIterator = cigar.getCigarElements().iterator();
            CigarElement cigarElement = cigarElementIterator.next();
            // Skip all leading hard clips
            while (cigarElement.getOperator() == CigarOperator.HARD_CLIP) {
                totalHardClipCount += cigarElement.getLength();
                if (cigarElementIterator.hasNext())
                    cigarElement = cigarElementIterator.next();
                else
                    throw new ReviewedStingException("Read is entirely hardclipped, shouldn't be trying to clip it's cigar string");
            }
            // keep clipping until we hit stop
            while (index <= stop) {
                int shift = 0;
                if (cigarElement.getOperator().consumesReadBases())
                    shift = cigarElement.getLength();

                // we're still clipping or just finished perfectly
                if (index + shift == stop + 1) {
                    newCigar.add(new CigarElement(totalHardClipCount, CigarOperator.HARD_CLIP));
                }
                // element goes beyond what we need to clip
                else if (index + shift > stop + 1) {
                    int elementLengthAfterChopping = cigarElement.getLength() - (stop - index + 1);
                    newCigar.add(new CigarElement(totalHardClipCount, CigarOperator.HARD_CLIP));
                    newCigar.add(new CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
                }
                index += shift;
                if (index <= stop && cigarElementIterator.hasNext())
                    cigarElement = cigarElementIterator.next();
            }

            // add the remaining cigar elements
            while (cigarElementIterator.hasNext()) {
                cigarElement = cigarElementIterator.next();
                newCigar.add(new CigarElement(cigarElement.getLength(), cigarElement.getOperator()));
            }
        }

        // hard clip the end of the cigar string
        else {
            Iterator<CigarElement> cigarElementIterator = cigar.getCigarElements().iterator();
            CigarElement cigarElement = cigarElementIterator.next();

            // Keep marching on until we find the start
            while (index < start) {
                int shift = 0;
                if (cigarElement.getOperator().consumesReadBases())
                    shift = cigarElement.getLength();

                // we haven't gotten to the start yet, keep everything as is.
                if (index + shift < start)
                    newCigar.add(new CigarElement(cigarElement.getLength(), cigarElement.getOperator()));

                // element goes beyond our clip starting position
                else {
                    int elementLengthAfterChopping = start - index;
                    // if this last element is a HARD CLIP operator, just merge it with our hard clip operator to be added later
                    if (cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                        totalHardClipCount += elementLengthAfterChopping;
                    // otherwise, maintain what's left of this last operator
                    else
                        newCigar.add(new CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
                }
                index += shift;
                if (index < start && cigarElementIterator.hasNext())
                    cigarElement = cigarElementIterator.next();
            }
            newCigar.add(new CigarElement(totalHardClipCount, CigarOperator.HARD_CLIP));
        }
        return newCigar;
    }

    private int calculateAlignmentStartShift(Cigar oldCigar, Cigar newCigar) {
        int shift = 0;

        // Rewind to previous start (by counting everything that was already clipped in this read)
        for (CigarElement cigarElement : oldCigar.getCigarElements()) {
            if (!cigarElement.getOperator().consumesReferenceBases())
                shift -= cigarElement.getLength();
            else
                break;
        }

        // Advance to new start (by counting everything new that has been clipped )
        for (CigarElement cigarElement : newCigar.getCigarElements()) {
            if (!cigarElement.getOperator().consumesReferenceBases())
                shift += cigarElement.getLength();
            else
                break;
        }
        return shift;
    }
}
