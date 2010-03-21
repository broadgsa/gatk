package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;

import java.util.Collection;

/**
 * Walks over all pairs/collections of reads in a BAM file sorted by
 * read name.
 *
 * @author mhanna
 * @version 0.1
 */
@Requires({DataSource.READS})
public abstract class ReadPairWalker<MapType,ReduceType> extends Walker<MapType,ReduceType> {

    /**
     * Optionally filters out read pairs.
     * @param reads collections of all reads with the same read name.
     * @return True to process the reads with map/reduce; false otherwise.
     */
    public boolean filter(Collection<SAMRecord> reads) {
        // Keep all pairs by default.
        return true;
    }

    /**
     * Maps a read pair to a given reduce of type MapType.  Semantics determined by subclasser.
     * @param reads Collection of reads having the same name.
     * @return Semantics defined by implementer.
     */
    public abstract MapType map(Collection<SAMRecord> reads);

    // Given result of map function
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);

}
