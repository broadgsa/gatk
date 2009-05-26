package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.List;

/**
 *
 * User: aaron
 * Date: May 14, 2009
 * Time: 3:28:50 PM
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
 * @date May 14, 2009
 * <p/>
 * Class LocusWindowShardStrategy
 * <p/>
 * This function knows how to shard on a genome loc boundry.  It guarantees
 * a one-to-one mapping between a GenomeLoc and shard. 
 */
public class LocusIntervalShardStrategy extends LocusShardStrategy {
    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic       the seq dictionary
     * @param intervals file
     */
    LocusIntervalShardStrategy(SAMSequenceDictionary dic, GenomeLocSortedSet intervals) {
        super(dic, intervals);
    }

    /**
     * This is how the various shards strategies implements their approach, adjusting this value
     *
     * @return the next shard size
     */
    protected long nextShardSize() {
        long nextSize = this.getCurrentInterval().getStop() - this.getCurrentInterval().getStart();
        return nextSize;
    }

    /**
     * set the next shards size
     *
     * @param size adjust the next size to this
     */
    public void adjustNextShardSize(long size) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

}
