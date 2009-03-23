package org.broadinstitute.sting.gatk;

import net.sf.samtools.SAMRecord;

import java.util.List;
import java.lang.ref.Reference;

import org.broadinstitute.sting.utils.GenomeLoc;
import edu.mit.broad.picard.reference.ReferenceSequence;

/**
 * Useful class for forwarding on locusContext data from this iterator
 * 
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:01:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class LocusContext {
    private GenomeLoc loc = null;
    private List<SAMRecord> reads = null;
    private List<Integer> offsets = null;
    private ReferenceSequence refContig = null;

    /**
     * Create a new LocusContext object
     *
     * @param loc
     * @param reads
     * @param offsets
     */
    public LocusContext(GenomeLoc loc, List<SAMRecord> reads, List<Integer> offsets) {
        this.loc = loc;
        this.reads = reads;
        this.offsets = offsets;
    }

    /**
     * get all of the reads within this context
     * 
     * @return
     */
    public List<SAMRecord> getReads() { return reads; }

    /**
     * Are there any reads associated with this locus?
     *
     * @return
     */
    public boolean hasReads() {
        return reads != null;
    }

    /**
     * How many reads cover this locus?
     * @return
     */
    public int numReads() {
        assert( reads != null );
        return reads.size();
    }

    /**
     * get a list of the equivalent positions within in the reads at Pos
     *
     * @return
     */
    public List<Integer> getOffsets() {
        return offsets;
    }

    public String getContig() { return getLocation().getContig(); }
    public long getPosition() { return getLocation().getStart(); }
    public GenomeLoc getLocation() { return loc; }

    /**
     * Returns the entire reference sequence contig associated with these reads
     *
     * @return ReferenceSequence object, or null if unavailable
     */
    public ReferenceSequence getReferenceContig() {
        return refContig;
    }

    /**
     * @return True if reference sequence contig is available
     */
    public boolean hasReferenceContig() {
        return refContig != null;
    }

    /**
     * Sets the reference sequence for this locus to contig
     * 
     * @param contig
     */
    public void setReferenceContig(final ReferenceSequence contig) {
        refContig = contig;
    }
}
