package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 7:09:22 PM
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
 * Class ShardStrategyFactory
 * <p/>
 * The Shard Strategy Factory,  use this class to create and transfer shard strategies
 * between different approaches.
 */
public class ShardStrategyFactory {
    public enum SHATTER_STRATEGY {
        MONOLITHIC,   // Put all of the available data into one shard.
        LOCUS_EXPERIMENTAL,
        READS_EXPERIMENTAL
    }

    /**
     * get a new shatter strategy
     *
     * @param readsDataSource     File pointer to BAM.
     * @param referenceDataSource File pointer to reference.
     * @param strat        what's our strategy - SHATTER_STRATEGY type
     * @param dic          the seq dictionary
     * @param startingSize the starting size
     * @return a shard strategy capable of dividing input data into shards.
     */
    static public ShardStrategy shatter(SAMDataSource readsDataSource, IndexedFastaSequenceFile referenceDataSource, SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize, GenomeLocParser genomeLocParser) {
        return ShardStrategyFactory.shatter(readsDataSource, referenceDataSource, strat, dic, startingSize, genomeLocParser, -1L);
    }

    /**
     * get a new shatter strategy
     *
     * @param readsDataSource     File pointer to BAM.
     * @param referenceDataSource File pointer to reference.
     * @param strat               what's our strategy - SHATTER_STRATEGY type
     * @param dic                 the seq dictionary
     * @param startingSize        the starting size
     * @return a shard strategy capable of dividing input data into shards.
     */
    static public ShardStrategy shatter(SAMDataSource readsDataSource, IndexedFastaSequenceFile referenceDataSource, SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize, GenomeLocParser genomeLocParser, long limitByCount) {
        switch (strat) {
            case LOCUS_EXPERIMENTAL:
                return new LocusShardStrategy(readsDataSource,referenceDataSource,genomeLocParser,null);
            case READS_EXPERIMENTAL:
                return new ReadShardStrategy(genomeLocParser,readsDataSource,null);
            default:
                throw new ReviewedStingException("Strategy: " + strat + " isn't implemented for this type of shatter request");
        }

    }


    /**
     * get a new shatter strategy
     *
     * @param readsDataSource     File pointer to BAM.
     * @param referenceDataSource File pointer to reference.
     * @param strat        what's our strategy - SHATTER_STRATEGY type
     * @param dic          the seq dictionary
     * @param startingSize the starting size
     * @return a shard strategy capable of dividing input data into shards.
     */
    static public ShardStrategy shatter(SAMDataSource readsDataSource, IndexedFastaSequenceFile referenceDataSource, SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize, GenomeLocParser genomeLocParser, GenomeLocSortedSet lst) {
        return ShardStrategyFactory.shatter(readsDataSource, referenceDataSource, strat, dic, startingSize, genomeLocParser, lst, -1l);

    }

    /**
     * get a new shatter strategy
     *
     * @param readsDataSource The reads used to shatter this file.
     * @param referenceDataSource The reference used to shatter this file.
     * @param strat        what's our strategy - SHATTER_STRATEGY type
     * @param dic          the seq dictionary
     * @param startingSize the starting size
     * @return A strategy for shattering this data.
     */
    static public ShardStrategy shatter(SAMDataSource readsDataSource, IndexedFastaSequenceFile referenceDataSource, SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize, GenomeLocParser genomeLocParser, GenomeLocSortedSet lst, long limitDataCount) {
        switch (strat) {
            case LOCUS_EXPERIMENTAL:
                return new LocusShardStrategy(readsDataSource,referenceDataSource,genomeLocParser,lst);
            case READS_EXPERIMENTAL:
                return new ReadShardStrategy(genomeLocParser, readsDataSource,lst);
            default:
                throw new ReviewedStingException("Strategy: " + strat + " isn't implemented");
        }

    }

}
