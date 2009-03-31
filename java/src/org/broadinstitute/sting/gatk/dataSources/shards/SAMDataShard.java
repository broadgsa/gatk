package org.broadinstitute.sting.gatk.dataSources.shards;

import edu.mit.broad.picard.sam.MergingSamRecordIterator;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: aaronmckenna
 * Date: Mar 29, 2009
 * Time: 8:47:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class SAMDataShard implements DataShard {

    // our iterator
    final private MergingSamRecordIterator iterator;

    // divide by reads or by loci
    private boolean byReads = true;

    // iterator bounds limiter
    private int lengthCount = 0;
    private final int limiter;

    public SAMDataShard(MergingSamRecordIterator iterator, int limiter) {
        this.iterator = iterator;
        this.limiter = limiter;
    }

    public SAMDataShard(MergingSamRecordIterator iterator) {
        this.iterator = iterator;
        limiter = Integer.MAX_VALUE;
    }


    public boolean hasNext() {
        return iterator.hasNext() && lengthCount > limiter;
    }

    public SAMRecord next() {
        ++lengthCount;
        return iterator.next();
    }

    public void remove() {
        iterator.remove();
    }
}
