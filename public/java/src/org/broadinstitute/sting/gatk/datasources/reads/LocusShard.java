package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.SAMFileSpan;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;
import java.util.Map;

/**
 * Handles locus shards of BAM information.
 * @author aaron
 * @version 1.0
 * @date Apr 7, 2009
 */
public class LocusShard extends Shard {
    /**
     * Create a new locus shard, divided by index.
     * @param intervals List of intervals to process.
     * @param fileSpans File spans associated with that interval.
     */
    public LocusShard(GenomeLocParser parser, SAMDataSource dataSource, List<GenomeLoc> intervals, Map<SAMReaderID,SAMFileSpan> fileSpans) {
        super(parser, ShardType.LOCUS, intervals, dataSource, fileSpans, false);
    }

    /**
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */
    @Override
    public String toString() {
        return Utils.join(";",getGenomeLocs());
    }
}
