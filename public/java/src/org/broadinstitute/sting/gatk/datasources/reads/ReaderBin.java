package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.Bin;
import net.sf.samtools.BrowseableBAMIndex;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Feb 2, 2011
 * Time: 4:36:40 PM
 * To change this template use File | Settings | File Templates.
 */
class ReaderBin {
    public final SAMReaderID id;
    public final BrowseableBAMIndex index;
    public final int referenceSequence;
    public final Bin bin;

    public ReaderBin(final SAMReaderID id, final BrowseableBAMIndex index, final int referenceSequence, final Bin bin) {
        this.id = id;
        this.index = index;
        this.referenceSequence = referenceSequence;
        this.bin = bin;
    }

    public int getStart() {
        return index.getFirstLocusInBin(bin);
    }

    public int getStop() {
        return index.getLastLocusInBin(bin);
    }
}
