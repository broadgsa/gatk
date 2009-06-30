package org.broadinstitute.sting.gatk;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.*;

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
        //assert loc != null;
        //assert loc.getContig() != null;
        //assert reads != null;
        //assert offsets != null;

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
    public void setLocation(GenomeLoc loc) {
        this.loc = loc.clone();
    }
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

    public void downsampleToCoverage(int coverage) {
        if ( numReads() <= coverage )
            return;

        // randomly choose numbers corresponding to positions in the reads list
        Random generator = new Random();
        TreeSet positions = new TreeSet();
        int i = 0;
        while ( i < coverage ) {
            if (positions.add(new Integer(generator.nextInt(reads.size()))))
                i++;
        }

        ArrayList<SAMRecord> downsampledReads = new ArrayList<SAMRecord>();
        ArrayList<Integer> downsampledOffsets = new ArrayList<Integer>();
        Iterator positionIter = positions.iterator();
        Iterator<SAMRecord> readsIter = reads.iterator();
        Iterator<Integer> offsetsIter = offsets.iterator();
        int currentRead = 0;
        while ( positionIter.hasNext() ) {
            int nextReadToKeep = (Integer)positionIter.next();

            // fast-forward to the right read
            while ( currentRead < nextReadToKeep ) {
                readsIter.next();
                offsetsIter.next();
                currentRead++;
            }

            downsampledReads.add(readsIter.next());
            downsampledOffsets.add(offsetsIter.next());
            currentRead++;
        }

        reads = downsampledReads;
        offsets = downsampledOffsets;
    }
}
