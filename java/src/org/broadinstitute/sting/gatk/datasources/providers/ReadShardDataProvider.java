package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Collection;

import net.sf.picard.reference.IndexedFastaSequenceFile;

/**
 * Present data sharded by read to a traversal engine.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadShardDataProvider extends ShardDataProvider {
    /**
     * The raw collection of reads.
     */
    private final StingSAMIterator reads;

    /**
     * Create a data provider for the shard given the reads and reference.
     * @param shard The chunk of data over which traversals happen.
     * @param reference A getter for a section of the reference.
     */
    public ReadShardDataProvider(Shard shard, StingSAMIterator reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        super(shard,reference,rods);
        this.reads = reads;
    }

    /**
     * Can this data source provide reads?
     * @return True if reads are available, false otherwise.
     */
    public boolean hasReads() {
        return reads != null;
    }    

    /**
     * Gets an iterator over all the reads bound by this shard.
     * @return An iterator over all reads in this shard.
     */
    public StingSAMIterator getReadIterator() {
        return reads;
    }

    @Override
    public void close() {
        super.close();
        
        if(reads != null)
            reads.close();
    }

}
