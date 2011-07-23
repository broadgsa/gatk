package org.broadinstitute.sting.utils.clipreads;

import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * A simple collection of the clipping operations to apply to a read along with its read
 */
public class ReadClipper {
    SAMRecord read;
    List<ClippingOp> ops = null;

    /**
     * We didn't do any clipping work on this read, just leave everything as a default
     *
     * @param read
     */
    public ReadClipper(final SAMRecord read) {
        this.read = read;
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
        return ops != null;
    }

    public SAMRecord getRead() {
        return read;
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
                return clippedRead;
            } catch (CloneNotSupportedException e) {
                throw new RuntimeException(e); // this should never happen
            }
        }
    }
}
