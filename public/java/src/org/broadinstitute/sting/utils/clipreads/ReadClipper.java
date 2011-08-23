package org.broadinstitute.sting.utils.clipreads;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import org.broad.tribble.util.PositionalStream;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.jets3t.service.multi.ThreadedStorageService;

import java.util.ArrayList;
import java.util.List;

/**
 * A simple collection of the clipping operations to apply to a read along with its read
 */
public class ReadClipper {
    SAMRecord read;
    boolean wasClipped;
    List<ClippingOp> ops = null;

    /**
     * We didn't do any clipping work on this read, just leave everything as a default
     *
     * @param read
     */
    public ReadClipper(final SAMRecord read) {
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

    public SAMRecord getRead() {
        return read;
    }

    public SAMRecord hardClipByReferenceCoordinates(int refStart, int refStop) {
        int start = (refStart < 0) ? 0 : ReadUtils.getReadCoordinateForReferenceCoordinate(read, refStart);
        int stop =  (refStop  < 0) ? read.getReadLength() - 1 : ReadUtils.getReadCoordinateForReferenceCoordinate(read, refStop);

        System.out.println("Clipping start/stop: " + start + "/" + stop);
        this.addOp(new ClippingOp(start, stop));
        return clipRead(ClippingRepresentation.HARDCLIP_BASES);
    }

    public SAMRecord hardClipByReadCoordinates(int start, int stop) {
        this.addOp(new ClippingOp(start, stop));
        return clipRead(ClippingRepresentation.HARDCLIP_BASES);
    }

    public SAMRecord hardClipBothEndsByReferenceCoordinates(int left, int right) {
        this.read = hardClipByReferenceCoordinates(-1, left);
        this.ops = null; // reset the operations
        return hardClipByReferenceCoordinates(right, -1);
    }

    public SAMRecord hardClipLowQualEnds(byte lowQual) {
        byte [] quals = read.getBaseQualities();
        int leftClipIndex = 0;
        int rightClipIndex = read.getReadLength() - 1;

        // check how far we can clip both sides
        while (rightClipIndex >= 0 && quals[rightClipIndex] <= lowQual) rightClipIndex--;
        while (leftClipIndex < read.getReadLength() && quals[leftClipIndex] <= lowQual) leftClipIndex++;

        // if the entire read should be clipped, then return an empty read. (--todo: maybe null is better? testing this for now)
        if (leftClipIndex > rightClipIndex)
            return (new SAMRecord(read.getHeader()));

        if (rightClipIndex < read.getReadLength() - 1) {
            this.addOp(new ClippingOp(rightClipIndex + 1, read.getReadLength() - 1));
        }
        if (leftClipIndex > 0 ) {
            this.addOp(new ClippingOp(0, leftClipIndex - 1));
        }
        return this.clipRead(ClippingRepresentation.HARDCLIP_BASES);
    }

    /**
     * Return a new read corresponding to this.read that's been clipped according to ops, if any are present.
     *
     * @param algorithm
     * @return
     */
    public SAMRecord clipRead(ClippingRepresentation algorithm) {
        if (ops == null)
            return getRead();
        else {
            try {
                SAMRecord clippedRead = (SAMRecord) read.clone();
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
}
