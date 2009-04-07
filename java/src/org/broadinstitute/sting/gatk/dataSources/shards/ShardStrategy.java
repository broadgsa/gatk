package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Iterator;
/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 11:23:17 AM
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
 * @date Apr 6, 2009
 * <p/>
 * Interface Shard
 * <p/>
 * The shard interface, which controls how data is divided
 */
public abstract class ShardStrategy implements Iterator<Shard>, Iterable<Shard> {

    // this stores the seq dictionary, which is a reference for the
    // lengths and names of contigs, which you need to generate an iterative stratagy
    protected final SAMSequenceDictionary dic;

    // the current genome location
    protected GenomeLoc mLoc = null;

    // current seq location
    protected int seqLoc = 0;

    // the actual last size; this can change based on contig endings
    protected long lastGenomeLocSize = 0;

    // do we have another contig?
    private boolean nextContig = false;


    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic the seq dictionary
     */
    ShardStrategy(SAMSequenceDictionary dic) {
        this.dic = dic;
        mLoc = new GenomeLoc(dic.getSequence(0).getSequenceName(), 0, 0);
        if (dic.getSequences().size() > 0) {
            nextContig = true;
        }
    }

    /**
     * the copy constructor,
     *
     * @param old the old strategy
     */
    ShardStrategy(ShardStrategy old) {
        this.dic = old.dic;
        this.mLoc = old.mLoc;
        this.seqLoc = old.seqLoc;
        this.lastGenomeLocSize = old.lastGenomeLocSize;
        this.nextContig = old.nextContig;
    }

    /**
     *
     * Abstract methods that each strategy has to implement
     *
     */

    /**
     * set the next shards size
     *
     * @param size adjust the next size to this
     */
    public abstract void adjustNextShardSize(long size);


    /**
     * This is how the various shards strategies implements their approach
     *
     * @return the next shard size
     */
    abstract long nextShardSize();


    /**
     *
     * Concrete methods that each strategy does not have to implement
     *
     */

    /**
     * get the next shard, based on the return size of nextShardSize
     *
     * @return
     */
    public Shard next() {
        // lets get some background info on the problem
        long length = dic.getSequence(seqLoc).getSequenceLength();
        long proposedSize = nextShardSize();
        long nextStart = mLoc.getStop() + 1;
        // can we fit it into the current seq size?
        if (nextStart + proposedSize < length) {
            lastGenomeLocSize = proposedSize;
            mLoc = new GenomeLoc(dic.getSequence(seqLoc).getSequenceName(), nextStart, nextStart + proposedSize);
            return Shard.toShard(new GenomeLoc(dic.getSequence(seqLoc).getSequenceName(), nextStart, nextStart + proposedSize));
        }
        // else we can't make it in the current location, we have to stitch one together
        else {
            lastGenomeLocSize = nextStart + proposedSize - length;


            // move to the next contig
            jumpContig();
            return Shard.toShard(new GenomeLoc(dic.getSequence(seqLoc).getSequenceName(), nextStart, lastGenomeLocSize));
        }

    }

    /** jump to the next contig */
    private void jumpContig() {
        ++seqLoc;
        if (dic.getSequences().size() <= seqLoc) {
            nextContig = false;
            return;
        }

        // the next sequence should start at the begining of the next contig
        mLoc = new GenomeLoc(dic.getSequence(seqLoc).getSequenceName(), 0, 0);

    }

    /**
     * is there another GenomeLoc to get?
     *
     * @return
     */
    public boolean hasNext() {
        return nextContig;
    }

    /** we don't support remove */
    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a shard iterator!");
    }


    /**
     * to be for-each(able), we must implement this method
     *
     * @return
     */
    public Iterator<Shard> iterator() {
        return this;
    }


}
