package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.GATKException;

import java.util.ArrayList;
import java.util.List;
import java.util.Collection;

import net.sf.picard.reference.IndexedFastaSequenceFile;
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
public abstract class ShardDataProvider {
    /**
     * An ArrayList of all the views that are examining this data.
     */
    private List<View> registeredViews = new ArrayList<View>();

    /**
     * The shard over which we're providing data.
     */
    private final Shard shard;

    /**
     * Provider of reference data for this particular shard.
     */
    private final IndexedFastaSequenceFile reference;

    /**
     * Sources of reference-ordered data.
     */
    private final Collection<ReferenceOrderedDataSource> referenceOrderedData;

    /**
     * Retrieves the shard associated with this data provider.
     * @return The shard associated with this data provider.
     */
    public Shard getShard() {
        return shard;
    }

    /**
     * Can this data source provide reference information?
     * @return True if possible, false otherwise.
     */
    public boolean hasReference() {
        return reference != null;
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
    Collection<ReferenceOrderedDataSource> getReferenceOrderedData() {
        return referenceOrderedData;        
    }

    /**
     * Create a data provider for the shard given the reads and reference.
     * @param shard The chunk of data over which traversals happen.
     * @param reference A getter for a section of the reference.
     */
    public ShardDataProvider(Shard shard,IndexedFastaSequenceFile reference,Collection<ReferenceOrderedDataSource> rods) {
        this.shard = shard;
        this.reference = reference;
        this.referenceOrderedData = rods;
    }

    /**
     * Skeletal, package protected constructor for unit tests which require a ShardDataProvider.
     * @param shard the shard
     */
    ShardDataProvider(Shard shard) {
        this(shard,null,null);
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
                    throw new GATKException(String.format("Tried to register two conflicting views: %s and %s",
                                                           registeredView.getClass().getSimpleName(),
                                                           view.getClass().getSimpleName()));
            }
        }

        // Check whether this class has any objection to any other classes.
        for( Class<? extends View> conflict: view.getConflictingViews() ) {
            for( View registeredView: registeredViews ) {
                if( conflict.isInstance(registeredView) )
                    throw new GATKException(String.format("Tried to register two conflicting views: %s and %s",
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

        // Explicitly purge registered views to ensure that we don't end up with circular references
        // to views, which can in turn hold state.
        registeredViews.clear();

        if(shard != null)
            shard.close();
    }

    @Override
    public String toString() {
        return shard.toString();
    }
}
