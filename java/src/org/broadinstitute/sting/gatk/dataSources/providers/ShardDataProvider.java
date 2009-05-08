package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ReadShard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLoc;
/**
 * User: hanna
 * Date: May 8, 2009
 * Time: 3:09:57 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * An umbrella class that examines the data passed to the microscheduler and
 * tries to assemble as much as possible with it. 
 */
public class ShardDataProvider {
    /**
     * The raw collection of reads.
     */
    private final StingSAMIterator reads;

    /**
     * Information about the locus.  Can be accessed mutually exclusively
     * with the reads iterator.
     */
    private final LocusContextProvider locusContextProvider;

    /**
     * Provider of reference data for this particular shard.
     */
    private final ReferenceProvider referenceProvider;

    /**
     * Users should only drive using the reads directly or using the locus.
     * Perhaps in the future we can support direct simultaneous access to both. 
     */
    private enum LastAccessType { READS, LOCUS };
    private LastAccessType lastAccessed = null;

    /**
     * Can this data source provide reads?
     * @return True if reads are available, false otherwise.
     */
    public boolean hasReads() {
        return reads != null;
    }

    /**
     * Can this data source provide a locus context?
     * @return True if possible, false otherwise.
     */
    public boolean hasLocusContext() {
        return locusContextProvider != null;
    }

    /**
     * Can this data source provide reference information?
     * @return True if possible, false otherwise.
     */
    public boolean hasReference() {
        return referenceProvider != null;
    }

    /**
     * Gets an iterator over all the reads bound by this shard.
     * WARNING: Right now, this cannot be concurrently accessed with getLocusContext().
     * @return An iterator over all reads in this shard.
     */
    public StingSAMIterator getReadIterator() {
        if( LastAccessType.LOCUS.equals(lastAccessed) )
            throw new UnsupportedOperationException("Cannot mix direct access to reads with access to locus context");
        lastAccessed = LastAccessType.READS;
        return reads;
    }

    /**
     * Gets a locus context for a particular region on the genome.
     * WARNING: Right now, this cannot be concurrently accessed with the read iterator.
     * WARNING: Right now, accesses must be sequential along the genome.  
     * @param genomeLoc The location for which to determine the context.
     * @return The context associated with this location.
     */
    public LocusContext getLocusContext( GenomeLoc genomeLoc ) {
        if( LastAccessType.READS.equals(lastAccessed) )
            throw new UnsupportedOperationException("Cannot mix direct access to reads with access to locus context");
        lastAccessed = LastAccessType.LOCUS;
        return locusContextProvider.getLocusContext( genomeLoc );
    }

    /**
     * Gets the reference base associated with this particular point on the genome.
     * @param genomeLoc Region for which to retrieve the base.  GenomeLoc must represent a 1-base region.
     * @return The base at the position represented by this genomeLoc.
     */
    public char getReferenceBase( GenomeLoc genomeLoc ) {
        return referenceProvider.getReferenceBase(genomeLoc);        
    }

    /**
     * Create a data provider for the shard given the reads and reference.
     * @param shard The chunk of data over which traversals happen.
     * @param reads A window into the reads for a given region.                                                
     * @param reference A getter for a section of the reference.
     */
    public ShardDataProvider( Shard shard, SAMDataSource reads, IndexedFastaSequenceFile reference ) {
        // Provide basic reads information.
        this.reads = reads.seek( shard );
        // Watch out!  the locus context provider will start prefetching data off the queue.  Only create this
        // if absolutely necessary.
        this.locusContextProvider = !(shard instanceof ReadShard) ? new LocusContextProvider(this.reads) : null;
        this.referenceProvider = (reference != null) ? new ReferenceProvider(reference,shard) : null;
    }

    /**
     * Retire this shard.
     */
    public void close() {
        reads.close();
    }
}
