package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.StingException;

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
 */
public abstract class LocusShardStrategy implements ShardStrategy {

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

    /** our interal list * */
    private GenomeLocSortedSet intervals = null;

    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(LocusShardStrategy.class);


    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic the seq dictionary
     */
    LocusShardStrategy(SAMSequenceDictionary dic) {
        this.dic = dic;
        mLoc = new GenomeLoc(0, 0, 0);
        if (dic.getSequences().size() > 0) {
            nextContig = true;
        }
    }

    /**
     * the copy constructor,
     *
     * @param old the old strategy
     */
    LocusShardStrategy(LocusShardStrategy old) {
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
    LocusShardStrategy(SAMSequenceDictionary dic, GenomeLocSortedSet intervals) {
        this.dic = dic;
        this.intervals = intervals.clone();
        // set the starting point to the beginning interval
        if (intervals.size() < 1) {
            throw new IllegalArgumentException("Interval files must contain at least one interval");
        }
        GenomeLoc loc = intervals.iterator().next();
        mLoc = new GenomeLoc(loc.getContig(), loc.getStart() - 1, loc.getStart() - 1);
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
     * This is how the various shards strategies implements their approach, adjusting this value
     *
     * @return the next shard size
     */
    protected abstract long nextShardSize();


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

        // if we don't have an interval set, use the non interval based approach.  Simple, eh?
        if (this.intervals == null) {
            return nonIntervaledNext(length, proposedSize, nextStart);
        } else {
            return intervaledNext(proposedSize);
        }

    }

    /**
     * Interval based next processing
     *
     * @param proposedSize the proposed size
     *
     * @return the shard that represents this data
     */
    private Shard intervaledNext(long proposedSize) {
        if ((this.intervals == null) || (intervals.isEmpty())) {
            throw new StingException("LocusShardStrategy: genomic regions list is empty in next() function.");
        }

        // get the first region in the list
        GenomeLoc loc = intervals.iterator().next();

        if (loc.getStop() - loc.getStart() <= proposedSize) {
            intervals.removeRegion(loc);
            return new IntervalReadShard(loc);
        } else {
            GenomeLoc subLoc = new GenomeLoc(loc.getContigIndex(), loc.getStart(), loc.getStart() + proposedSize - 1);
            intervals.removeRegion(subLoc);
            return new IntervalReadShard(subLoc);
        }
    }

    /**
     * Get the next shard, if we don't have intervals to traverse over
     *
     * @param length       the length of the contig
     * @param proposedSize the proposed size
     * @param nextStart    the next start location
     *
     * @return the shard to return to the user
     */
    private Shard nonIntervaledNext(long length, long proposedSize, long nextStart) {
        // can we fit it into the current seq size?
        if (nextStart + proposedSize - 1 < length) {
            lastGenomeLocSize = proposedSize;
            mLoc = new GenomeLoc(dic.getSequence(seqLoc).getSequenceIndex(), nextStart, nextStart + proposedSize - 1);
            return LocusShard.toShard(new GenomeLoc(dic.getSequence(seqLoc).getSequenceIndex(), nextStart, nextStart + proposedSize - 1));
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
            Shard ret = LocusShard.toShard(new GenomeLoc(dic.getSequence(seqLoc).getSequenceIndex(), nextStart, nextStart + lastGenomeLocSize - 1));

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
        logger.debug("Next contig, index = " + dic.getSequence(seqLoc).getSequenceIndex());
        mLoc = new GenomeLoc(dic.getSequence(seqLoc).getSequenceIndex(), 0, 0);


    }

    /**
     * is there another GenomeLoc to get?
     *
     * @return
     */
    public boolean hasNext() {
        // if we don't have an interval file, use the non interval based approach.
        if (this.intervals == null) {
            return nextContig;
        } else {
            return (this.intervals.size() > 0);
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

    /**
     * this allows a shard strategy to get the current interval.  It's kind of a hack, but for the
     * locusWindowShardStrategy it was the best approach.
     *
     * @return
     */
    protected GenomeLoc getCurrentInterval() {
        if (this.intervals == null || intervals.size() < 1) {
            return null;
        }
        return intervals.iterator().next();
    }

}
