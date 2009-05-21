package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;
import java.util.ArrayList;
/**
 * User: hanna
 * Date: May 21, 2009
 * Time: 2:49:17 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A view into the reference-ordered data in the provider.
 */
public class ReferenceOrderedView implements View {
    /**
     * The provider that's supplying our backing data.
     */
    private final ShardDataProvider provider;

    /**
     * The data sources along with their current states.
     */
    private List<ReferenceOrderedDataState> states = new ArrayList<ReferenceOrderedDataState>();

    /**
     * Create a new view of reference-ordered data.
     * @param provider
     */
    public ReferenceOrderedView( ShardDataProvider provider ) {
        this.provider = provider;
        for( ReferenceOrderedDataSource dataSource: provider.getReferenceOrderedData() )
            states.add( new ReferenceOrderedDataState( dataSource, (ReferenceOrderedData.RODIterator)dataSource.seek(provider.getShard()) ) );

        provider.register(this);        
    }

    /**
     * Gets an object which can track the reference-ordered data at every locus.
     * @param loc Locus at which to track.
     * @return A tracker containing information about this locus.
     */
    public RefMetaDataTracker getReferenceOrderedDataAtLocus( GenomeLoc loc ) {
        RefMetaDataTracker tracks = new RefMetaDataTracker();
        for (ReferenceOrderedDataState state: states )
            tracks.bind( state.dataSource.getName(), state.iterator.seekForward(loc) );
        return tracks;
    }

    /**
     * Closes the current view.
     */
    public void close() {
        for( ReferenceOrderedDataState state: states )
            state.dataSource.close( state.iterator );

        // Clear out the existing data so that post-close() accesses to this data will fail-fast.
        states = null;
    }

    /**
     * Models the traversal state of a given ROD lane.
     */
    private class ReferenceOrderedDataState {
        public final ReferenceOrderedDataSource dataSource;
        public final ReferenceOrderedData.RODIterator iterator;

        public ReferenceOrderedDataState( ReferenceOrderedDataSource dataSource, ReferenceOrderedData.RODIterator iterator ) {
            this.dataSource = dataSource;
            this.iterator = iterator;
        }
    }
}
