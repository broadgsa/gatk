package org.broadinstitute.sting.gatk.datasources.providers;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Collection;

/**
 * Presents data sharded by locus to the traversal engine.
 *
 * @author mhanna
 * @version 0.1
 */
public class LocusShardDataProvider extends ShardDataProvider {
    /**
     * Information about the source of the read data.
     */
    private final ReadProperties sourceInfo;

    /**
     * The parser, used to create and build new GenomeLocs.
     */
    private final GenomeLocParser genomeLocParser;

    /**
     * The particular locus for which data is provided.  Should be contained within shard.getGenomeLocs().
     */
    private final GenomeLoc locus;

    /**
     * The raw collection of reads.
     */
    private final LocusIterator locusIterator;

    /**
     * Create a data provider for the shard given the reads and reference.
     * @param shard The chunk of data over which traversals happen.
     * @param reference A getter for a section of the reference.
     */
    public LocusShardDataProvider(Shard shard, ReadProperties sourceInfo, GenomeLocParser genomeLocParser, GenomeLoc locus, LocusIterator locusIterator, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        super(shard,genomeLocParser,reference,rods);
        this.sourceInfo = sourceInfo;
        this.genomeLocParser = genomeLocParser;
        this.locus = locus;
        this.locusIterator = locusIterator;
    }

    /**
     * Returns information about the source of the reads.
     * @return Info about the source of the reads.
     */
    public ReadProperties getSourceInfo() {
        return sourceInfo;
    }

    /**
     * Gets the locus associated with this shard data provider.
     * @return The locus.
     */
    public GenomeLoc getLocus() {
        return locus;
    }

    /**
     * Gets an iterator over all the reads bound by this shard.
     * @return An iterator over all reads in this shard.
     */
    public LocusIterator getLocusIterator() {
        return locusIterator;
    }        

    @Override
    public void close() {
        super.close();
    }
}
