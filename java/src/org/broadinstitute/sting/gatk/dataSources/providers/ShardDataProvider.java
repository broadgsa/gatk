package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.samtools.SAMRecord;
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
     * The shard over which we're providing data.
     */
    private final Shard shard;

    /**
     * The raw collection of reads.
     */
    private final StingSAMIterator reads;

    /**
     * Provider of reference data for this particular shard.
     */
    private final ReferenceProvider referenceProvider;

    /**
     * Retrieves the shard associated with this data provider.
     * @return The shard associated with this data provider.
     */
    public Shard getShard() {
        return shard;
    }

    /**
     * Can this data source provide reads?
     * @return True if reads are available, false otherwise.
     */
    public boolean hasReads() {
        return reads != null;
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
        return reads;
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
     * Gets the reference sequence, as a char[], for the provided read.
     * @param read the read to fetch the reference sequence for
     * @return a char string of bases representing the reference sequence mapped to passed in read
     */
    public char[] getReferenceForRead( SAMRecord read ) {
        return referenceProvider.getReferenceBases(read);
    }

    /**
     * Create a data provider for the shard given the reads and reference.
     * @param shard The chunk of data over which traversals happen.
     * @param reads A window into the reads for a given region.                                                
     * @param reference A getter for a section of the reference.
     */
    public ShardDataProvider( Shard shard, SAMDataSource reads, IndexedFastaSequenceFile reference ) {
        this.shard = shard;
        // Provide basic reads information.
        this.reads = reads.seek( shard );
        this.referenceProvider = (reference != null) ? new ReferenceProvider(reference,shard) : null;
    }

    /**
     * Skeletal, package protected constructor for unit tests which require a ShardDataProvider.
     * @param shard the shard
     * @param reads reads iterator.
     */
    ShardDataProvider( Shard shard, StingSAMIterator reads ) {
        this.shard = shard;
        this.reads = reads;
        this.referenceProvider = null;
    }

    /**
     * Retire this shard.
     */
    public void close() {
        reads.close();
    }
}
