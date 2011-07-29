package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
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
public class ManagingReferenceOrderedView implements ReferenceOrderedView {
    /**
     * The data sources along with their current states.
     */
    private List<ReferenceOrderedDataState> states = new ArrayList<ReferenceOrderedDataState>();

    /**
     * Create a new view of reference-ordered data.
     * @param provider
     */
    public ManagingReferenceOrderedView( LocusShardDataProvider provider ) {
        for( ReferenceOrderedDataSource dataSource: provider.getReferenceOrderedData() )
            states.add(new ReferenceOrderedDataState(dataSource, dataSource.seek(provider.getLocus())));

        provider.register(this);
    }

    public Collection<Class<? extends View>> getConflictingViews() { return Collections.emptyList(); }

    /**
     * Gets an object which can track the reference-ordered data at every locus.
     * @param loc Locus at which to track.
     * @return A tracker containing information about this locus.
     */
    public RefMetaDataTracker getReferenceOrderedDataAtLocus( GenomeLoc loc, ReferenceContext referenceContext ) {
        List<RODRecordList> bindings = states.isEmpty() ? Collections.<RODRecordList>emptyList() : new ArrayList<RODRecordList>(states.size());

        for ( ReferenceOrderedDataState state: states )
            // todo -- warning, I removed the reference to the name from states
            bindings.add( state.iterator.seekForward(loc) );

        return new RefMetaDataTracker(bindings, referenceContext);
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
}

/**
 * Models the traversal state of a given ROD lane.
 */
class ReferenceOrderedDataState {
    public final ReferenceOrderedDataSource dataSource;
    public final LocationAwareSeekableRODIterator iterator;

    public ReferenceOrderedDataState( ReferenceOrderedDataSource dataSource, LocationAwareSeekableRODIterator iterator ) {
        this.dataSource = dataSource;
        this.iterator = iterator;
    }
}