package org.broadinstitute.sting.utils.clipreads;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.Vector;

/**
 * Represents a clip on a read.  It has a type (see the enum) along with a start and stop in the bases
 * of the read, plus an option extraInfo (useful for carrying info where needed).
 * <p/>
 * Also holds the critical apply function that actually execute the clipping operation on a provided read,
 * according to the wishes of the supplid ClippingAlgorithm enum.
 */
public class ClippingOp {
    public final ClippingType type;
    public final int start, stop; // inclusive
    public final Object extraInfo;

    public ClippingOp(int start, int stop) {
        this(null, start, stop, null);
    }

    public ClippingOp(ClippingType type, int start, int stop, Object extraInfo) {
        // todo -- remove type and extra info
        this.type = type;
        this.start = start;
        this.stop = stop;
        this.extraInfo = extraInfo;
    }

    public int getLength() {
        return stop - start + 1;
    }

    /**
     * Clips the bases in clippedRead according to this operation's start and stop.  Uses the clipping
     * representation used is the one provided by algorithm argument.
     *
     * @param algorithm
     * @param clippedRead
     */
    public SAMRecord apply(ClippingRepresentation algorithm, SAMRecord clippedRead) {
        //clippedRead.setReferenceIndex(1);
        byte[] quals = clippedRead.getBaseQualities();
        byte[] bases = clippedRead.getReadBases();

        switch (algorithm) {
            // important note:
            //   it's not safe to call read.getReadBases()[i] = 'N' or read.getBaseQualities()[i] = 0
            //   because you're not guaranteed to get a pointer to the actual array of bytes in the SAMRecord
            case WRITE_NS:
                for (int i = start; i <= stop; i++)
                    bases[i] = 'N';
                clippedRead.setReadBases(bases);
                break;
            case WRITE_Q0S:
                for (int i = start; i <= stop; i++)
                    quals[i] = 0;
                clippedRead.setBaseQualities(quals);
                break;
            case WRITE_NS_Q0S:
                for (int i = start; i <= stop; i++) {
                    bases[i] = 'N';
                    quals[i] = 0;
                }
                clippedRead.setReadBases(bases);
                clippedRead.setBaseQualities(quals);
                break;
            case HARDCLIP_BASES:
            case SOFTCLIP_BASES:
                if ( clippedRead.getReadUnmappedFlag() ) {
                    // we can't process unmapped reads
                    throw new UserException("Read Clipper cannot soft/hard clip unmapped reads");
                }

                //System.out.printf("%d %d %d%n", stop, start, clippedRead.getReadLength());
                int myStop = stop;
                if ( (stop + 1 - start) == clippedRead.getReadLength() ) {
                    // BAM representation issue -- we can't SOFTCLIP away all bases in a read, just leave it alone
                    //Walker.logger.info(String.format("Warning, read %s has all bases clip but this can't be represented with SOFTCLIP_BASES, just leaving it alone", clippedRead.getReadName()));
                    //break;
                    myStop--; // just decrement stop
                }

                if ( start > 0 && myStop != clippedRead.getReadLength() - 1 )
                    throw new RuntimeException(String.format("Cannot apply soft clipping operator to the middle of a read: %s to be clipped at %d-%d",
                            clippedRead.getReadName(), start, myStop));

                Cigar oldCigar = clippedRead.getCigar();

                int scLeft = 0, scRight = clippedRead.getReadLength();
                if ( start == 0 )
                    scLeft = myStop + 1;
                else
                    scRight = start;

                Cigar newCigar = softClip(oldCigar, scLeft, scRight);
                clippedRead.setCigar(newCigar);

                int newClippedStart = getNewAlignmentStartOffset(newCigar, oldCigar);
                int newStart = clippedRead.getAlignmentStart() + newClippedStart;
                clippedRead.setAlignmentStart(newStart);

                if ( algorithm == ClippingRepresentation.HARDCLIP_BASES )
                    clippedRead = ReadUtils.hardClipSoftClippedBases(clippedRead);
                //System.out.printf("%s clipping at %d %d / %d %d => %s and %d%n", oldCigar.toString(), start, stop, scLeft, scRight, newCigar.toString(), newStart);

                break;

            default:
                throw new IllegalStateException("Unexpected Clipping operator type " + algorithm);
        }

        return clippedRead;
    }

    /**
     * What is the type of a ClippingOp?
     */
    public enum ClippingType {
        LOW_Q_SCORES,
        WITHIN_CLIP_RANGE,
        MATCHES_CLIP_SEQ
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
                if (curLength > __startClipEnd && curLength < __endClipBegin) {
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
}
