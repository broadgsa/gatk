package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.Bin;

import java.util.HashMap;
import java.util.Map;

/**
 * Models a bin at which all BAM files in the merged input stream overlap.
 */
class BAMOverlap {
    public final int start;
    public final int stop;

    private final Map<SAMReaderID,Bin> bins = new HashMap<SAMReaderID,Bin>();

    public BAMOverlap(final int start, final int stop) {
        this.start = start;
        this.stop = stop;
    }

    public void addBin(final SAMReaderID id, final Bin bin) {
        bins.put(id,bin);
    }

    public Bin getBin(final SAMReaderID id) {
        return bins.get(id);
    }
}
