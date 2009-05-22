package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.NullSAMIterator;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.StingException;

import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
import java.io.File;
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
     * An ArrayList of all the views that are examining this data.
     */
    private List<View> registeredViews = new ArrayList<View>();

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
    private final IndexedFastaSequenceFile reference;

    /**
     * Sources of reference-ordered data.
     */
    private final List<ReferenceOrderedDataSource> referenceOrderedData;

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
        return reference != null;
    }

    /**
     * Gets an iterator over all the reads bound by this shard.
     * @return An iterator over all reads in this shard.
     */
    StingSAMIterator getReadIterator() {
        return reads;
    }

    /**
     * Gets a pointer into the given indexed fasta sequence file.
     * @return The indexed fasta sequence file.
     */
    IndexedFastaSequenceFile getReference() {
        return reference;        
    }

    /**
     * Gets a window into the reference-ordered data.  Package protected so that only
     * views can access it.
     * @return List of reference-ordered data sources.
     */
    List<ReferenceOrderedDataSource> getReferenceOrderedData() {
        return referenceOrderedData;        
    }

    /**
     * Create a data provider for the shard given the reads and reference.
     * @param shard The chunk of data over which traversals happen.
     * @param reads A window into the reads for a given region.
     * @param reference A getter for a section of the reference.
     */
    public ShardDataProvider( Shard shard, SAMDataSource reads, IndexedFastaSequenceFile reference, List<ReferenceOrderedDataSource> rods) {
        this.shard = shard;
        // Provide basic reads information.
        this.reads = (reads != null) ? reads.seek( shard ) : new NullSAMIterator(new Reads(new ArrayList<File>()));
        this.reference = reference;
        this.referenceOrderedData = rods;
    }

    /**
     * Skeletal, package protected constructor for unit tests which require a ShardDataProvider.
     * @param shard the shard
     * @param reads reads iterator.
     */
    ShardDataProvider( Shard shard, StingSAMIterator reads ) {
        this.shard = shard;
        this.reads = reads;
        this.reference = null;
        this.referenceOrderedData = null;
    }

    /**
     * Register this view with the shard provider, and make sure it has no conflicts with any other views.
     * @param view The new view.
     */
    void register( View view ) {
        // Check all registered classes to see whether a conflict exists.
        for( View registeredView: registeredViews ) {
            Collection<Class<? extends View>> conflicts = registeredView.getConflictingViews();
            for( Class<? extends View> conflict: conflicts ) {
                if( conflict.isInstance(view) )
                    throw new StingException(String.format("Tried to registered two conflicting views: %s and %s",
                                                           registeredView.getClass().getSimpleName(),
                                                           view.getClass().getSimpleName()));
            }
        }

        // Check whether this class has any objection to any other classes.
        for( Class<? extends View> conflict: view.getConflictingViews() ) {
            for( View registeredView: registeredViews ) {
                if( conflict.isInstance(registeredView) )
                    throw new StingException(String.format("Tried to registered two conflicting views: %s and %s",
                                                           registeredView.getClass().getSimpleName(),
                                                           view.getClass().getSimpleName()));
            }
        }

        this.registeredViews.add(view);
    }

    /**
     * Retire this shard.
     */
    public void close() {
        for( View view: registeredViews )
            view.close();
        reads.close();
    }
}
