package org.broadinstitute.sting.gatk.dataSources.shards;


import net.sf.samtools.SAMSequenceDictionary;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Iterator;
import java.util.List;
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

    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(ShardStrategy.class);

    /** our interal list * */
    private List<GenomeLoc> intervals = null;
    /** our interal list * */
    private int currentInterval = -1;

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
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic       the seq dictionary
     * @param intervals file
     */
    ShardStrategy(SAMSequenceDictionary dic, List<GenomeLoc> intervals) {
        this.dic = dic;
        this.intervals = intervals;
        this.currentInterval = 0;

        mLoc = new GenomeLoc(intervals.get(0).getContig(), intervals.get(0).getStart() - 1, intervals.get(0).getStart() - 1);
        if (dic.getSequences().size() > 0) {
            nextContig = true;
        }
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
     * @return the next shard
     */
    public Shard next() {

        // lets get some background info on the problem
        long length = dic.getSequence(seqLoc).getSequenceLength();
        long proposedSize = nextShardSize();
        long nextStart = mLoc.getStop() + 1;

        // if we don't have an interval file, use the non interval based approach.  Simple, eh?
        if (this.intervals == null) {
            return nonIntervaledNext(length, proposedSize, nextStart);
        } else {
            return intervaledNext(length, proposedSize, nextStart);
        }

    }

    private Shard intervaledNext(long length, long proposedSize, long nextStart) {
        // get the current genome location
        GenomeLoc loc = intervals.get(currentInterval);
        if (nextStart + proposedSize > loc.getStop()) {
            // we need to move the next interval
            proposedSize = loc.getStop() - nextStart;
            lastGenomeLocSize = proposedSize;

            // the next sequence should start at the begining of the next contig
            Shard ret = Shard.toShard(new GenomeLoc(intervals.get(currentInterval).getContig(), nextStart, nextStart + proposedSize - 1));

            ++currentInterval;
            if (intervals.size() > currentInterval) {
                mLoc = new GenomeLoc(intervals.get(currentInterval).getContig(), intervals.get(currentInterval).getStart() - 1, intervals.get(currentInterval).getStart() - 1);
            }
            return ret;// return 

        } else {
            // we need to move the next interval
            lastGenomeLocSize = proposedSize;

            // the next sequence should start at the begining of the next contig
            Shard ret = Shard.toShard(new GenomeLoc(intervals.get(currentInterval).getContig(), nextStart, nextStart + proposedSize - 1));

            mLoc = new GenomeLoc(intervals.get(currentInterval).getContig(), nextStart, nextStart + proposedSize - 1);

            return ret;// return
        }
    }

    /**
     * Get the next shard, if we don't have intervals to traverse over
     *
     * @param length       the length of the contig
     * @param proposedSize the proposed size
     * @param nextStart    the next start location
     * @return the shard to return to the user
     */
    private Shard nonIntervaledNext(long length, long proposedSize, long nextStart) {
        // can we fit it into the current seq size?
        if (nextStart + proposedSize - 1 < length) {
            lastGenomeLocSize = proposedSize;
            mLoc = new GenomeLoc(dic.getSequence(seqLoc).getSequenceName(), nextStart, nextStart + proposedSize - 1);
            return Shard.toShard(new GenomeLoc(dic.getSequence(seqLoc).getSequenceName(), nextStart, nextStart + proposedSize - 1));
        }
        // else we can't make it in the current location, we have to stitch one together
        else {
            // lets find out the remaining size of the current contig
            long overflow = nextStart + proposedSize - 1 - length;
            logger.debug("Overflow = " + overflow + " length: " + length);

            // set our last size counter to the remaining size
            lastGenomeLocSize = proposedSize - overflow;

            // move to the next contig
            // the next sequence should start at the begining of the next contig
            Shard ret = Shard.toShard(new GenomeLoc(dic.getSequence(seqLoc).getSequenceName(), nextStart, nextStart + lastGenomeLocSize));

            // now  jump ahead to the next contig
            jumpContig();

            // return the shard
            return ret;
        }
    }

    /** jump to the next contig */
    private void jumpContig() {
        ++seqLoc;

        if (!(seqLoc < dic.getSequences().size())) {
            nextContig = false;
            return;
        }
        logger.debug("Next contig, name = " + dic.getSequence(seqLoc).getSequenceName());
        mLoc = new GenomeLoc(dic.getSequence(seqLoc).getSequenceName(), 0, 0);


    }

    /**
     * is there another GenomeLoc to get?
     *
     * @return
     */
    public boolean hasNext() {
        // if we don't have an interval file, use the non interval based approach.  Simple, eh?
        if (this.intervals == null) {
            return nextContig;
        } else {
            return (this.currentInterval < this.intervals.size());
        }
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
