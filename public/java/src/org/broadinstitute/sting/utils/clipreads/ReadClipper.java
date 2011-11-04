package org.broadinstitute.sting.utils.clipreads;

import com.google.java.contract.Requires;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A simple collection of the clipping operations to apply to a read along with its read
 */
public class ReadClipper {
    GATKSAMRecord read;
    boolean wasClipped;
    List<ClippingOp> ops = null;

    /**
     * We didn't do any clipping work on this read, just leave everything as a default
     *
     * @param read
     */
    public ReadClipper(final GATKSAMRecord read) {
        this.read = read;
        this.wasClipped = false;
    }

    /**
     * Add another clipping operation to apply to this read
     *
     * @param op
     */
    public void addOp(ClippingOp op) {
        if (ops == null) ops = new ArrayList<ClippingOp>();
        ops.add(op);
    }

    public List<ClippingOp> getOps() {
        return ops;
    }

    public boolean wasClipped() {
        return wasClipped;
    }

    public GATKSAMRecord getRead() {
        return read;
    }

    public GATKSAMRecord hardClipByReferenceCoordinatesLeftTail(int refStop) {
        return hardClipByReferenceCoordinates(-1, refStop);
    }

    public GATKSAMRecord hardClipByReferenceCoordinatesRightTail(int refStart) {
        return hardClipByReferenceCoordinates(refStart, -1);
    }

    private int numDeletions(GATKSAMRecord read) {
        int result = 0;
        for (CigarElement e: read.getCigar().getCigarElements()) {
            if ( e.getOperator() == CigarOperator.DELETION || e.getOperator() == CigarOperator.D )
                result =+ e.getLength();
        }
        return result;
    }

    protected GATKSAMRecord hardClipByReferenceCoordinates(int refStart, int refStop) {
        int start = (refStart < 0) ? 0 : ReadUtils.getReadCoordinateForReferenceCoordinate(read, refStart, ReadUtils.ClippingTail.RIGHT_TAIL);
        int stop =  (refStop  < 0) ? read.getReadLength() - 1 : ReadUtils.getReadCoordinateForReferenceCoordinate(read, refStop, ReadUtils.ClippingTail.LEFT_TAIL);

        if (start < 0 || stop > read.getReadLength() - 1)
            throw new ReviewedStingException("Trying to clip before the start or after the end of a read");

        if ( start > stop )
            throw new ReviewedStingException("START > STOP -- this should never happen -- call Mauricio!");

        this.addOp(new ClippingOp(start, stop));
        GATKSAMRecord clippedRead = clipRead(ClippingRepresentation.HARDCLIP_BASES);
        this.ops = null;
        return clippedRead;
    }

    public GATKSAMRecord hardClipByReadCoordinates(int start, int stop) {
        this.addOp(new ClippingOp(start, stop));
        return clipRead(ClippingRepresentation.HARDCLIP_BASES);
    }

    @Requires("left <= right")
    public GATKSAMRecord hardClipBothEndsByReferenceCoordinates(int left, int right) {
        if (left == right)
            return new GATKSAMRecord(read.getHeader());
        GATKSAMRecord leftTailRead = hardClipByReferenceCoordinates(right, -1);

        // after clipping one tail, it is possible that the consequent hard clipping of adjacent deletions
        // make the left cut index no longer part of the read. In that case, clip the read entirely.
        if (left > leftTailRead.getAlignmentEnd())
            return new GATKSAMRecord(read.getHeader());

        ReadClipper clipper = new ReadClipper(leftTailRead);
        return clipper.hardClipByReferenceCoordinatesLeftTail(left);
    }

    public GATKSAMRecord hardClipLowQualEnds(byte lowQual) {
        byte [] quals = read.getBaseQualities();
        int leftClipIndex = 0;
        int rightClipIndex = read.getReadLength() - 1;

        // check how far we can clip both sides
        while (rightClipIndex >= 0 && quals[rightClipIndex] <= lowQual) rightClipIndex--;
        while (leftClipIndex < read.getReadLength() && quals[leftClipIndex] <= lowQual) leftClipIndex++;

        // if the entire read should be clipped, then return an empty read. (--todo: maybe null is better? testing this for now)
        if (leftClipIndex > rightClipIndex)
            return (new GATKSAMRecord(read.getHeader()));

        if (rightClipIndex < read.getReadLength() - 1) {
            this.addOp(new ClippingOp(rightClipIndex + 1, read.getReadLength() - 1));
        }
        if (leftClipIndex > 0 ) {
            this.addOp(new ClippingOp(0, leftClipIndex - 1));
        }
        return this.clipRead(ClippingRepresentation.HARDCLIP_BASES);
    }

    public GATKSAMRecord hardClipSoftClippedBases () {
        int readIndex = 0;
        int cutLeft = -1;            // first position to hard clip (inclusive)
        int cutRight = -1;           // first position to hard clip (inclusive)
        boolean rightTail = false;   // trigger to stop clipping the left tail and start cutting the right tail

        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() == CigarOperator.SOFT_CLIP) {
                if (rightTail) {
                    cutRight = readIndex;
                }
                else {
                    cutLeft = readIndex + cigarElement.getLength() - 1;
                }
            }
            else if (cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                rightTail = true;

            if (cigarElement.getOperator().consumesReadBases())
                readIndex += cigarElement.getLength();
        }

        // It is extremely important that we cut the end first otherwise the read coordinates change.
        if (cutRight >= 0)
            this.addOp(new ClippingOp(cutRight, read.getReadLength() - 1));
        if (cutLeft >= 0)
            this.addOp(new ClippingOp(0, cutLeft));

        return clipRead(ClippingRepresentation.HARDCLIP_BASES);
    }



    /**
     * Return a new read corresponding to this.read that's been clipped according to ops, if any are present.
     *
     * @param algorithm
     * @return
     */
    public GATKSAMRecord clipRead(ClippingRepresentation algorithm) {
        if (ops == null)
            return getRead();
        else {
            try {
                GATKSAMRecord clippedRead = (GATKSAMRecord) read.clone();
                for (ClippingOp op : getOps()) {
                    clippedRead = op.apply(algorithm, clippedRead);
                }
                wasClipped = true;
                return clippedRead;
            } catch (CloneNotSupportedException e) {
                throw new RuntimeException(e); // this should never happen
            }
        }
    }

    public GATKSAMRecord hardClipLeadingInsertions() {
        for(CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() != CigarOperator.HARD_CLIP && cigarElement.getOperator() != CigarOperator.SOFT_CLIP &&
                cigarElement.getOperator() != CigarOperator.INSERTION && cigarElement.getOperator() != CigarOperator.DELETION)
                break;

            else if (cigarElement.getOperator() == CigarOperator.INSERTION) {
                this.addOp(new ClippingOp(0, cigarElement.getLength() - 1));
            }

            else if (cigarElement.getOperator() == CigarOperator.DELETION) {
                throw new ReviewedStingException("No read should start with a deletion. Aligner bug?");
            }
        }
        return clipRead(ClippingRepresentation.HARDCLIP_BASES);
    }
}
