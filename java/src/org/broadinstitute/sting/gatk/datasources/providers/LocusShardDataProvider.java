package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.Reads;

import java.util.Collection;

import net.sf.picard.reference.IndexedFastaSequenceFile;

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
    private final Reads sourceInfo;

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
    public LocusShardDataProvider(Shard shard, Reads sourceInfo, GenomeLoc locus, LocusIterator locusIterator, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        super(shard,reference,rods);
        this.sourceInfo = sourceInfo;
        this.locus = locus;
        this.locusIterator = locusIterator;
    }

    /**
     * Returns information about the source of the reads.
     * @return Info about the source of the reads.
     */
    public Reads getSourceInfo() {
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
