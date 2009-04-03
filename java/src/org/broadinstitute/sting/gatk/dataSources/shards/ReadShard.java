package org.broadinstitute.sting.gatk.dataSources.shards;

import edu.mit.broad.picard.sam.MergingSamRecordIterator;
import org.broadinstitute.sting.gatk.dataSources.datum.ReadDatum;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.Arrays;

/**
 *
 * User: aaron
 * Date: Mar 30, 2009
 * Time: 5:45:51 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Mar 30, 2009
 * <p/>
 * Class ReadShard
 * <p/>
 * A read data shard.
 */
public class ReadShard implements DataShard {

    private MergingSamRecordIterator iterator;

    /**
     * create the data chunk with an iterator, and a limiter
     *
     * @param samIterator
     */
    public ReadShard(MergingSamRecordIterator samIterator) {
        this.iterator = samIterator;
    }

    /**
     * do we have a next data point
     *
     * @return true if we have a data point
     */
    public boolean hasNext() {
        return iterator.hasNext();
    }

    public ReadDatum next() {
        // get the read
        final SAMRecord read = iterator.next();

        // put the read into a list
        final List<SAMRecord> reads = Arrays.asList(read);

        // put together the genome location
        final GenomeLoc loc = GenomeLoc.genomicLocationOf(read);

        // Offset of a single read is always 0
        List<Integer> offsets = Arrays.asList(0);

        // create the locus
        final LocusContext locus = new LocusContext(loc, reads, offsets);

        // return the read datum
        return new ReadDatum(read, locus);
    }

    /** remove the current pointed to data source */
    public void remove() {
        iterator.remove();
    }
}
