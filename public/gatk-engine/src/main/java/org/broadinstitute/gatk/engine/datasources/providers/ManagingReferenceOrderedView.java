/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.datasources.providers;

import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;

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
    @Override
    public RefMetaDataTracker getReferenceOrderedDataAtLocus( GenomeLoc loc ) {
        if ( states.isEmpty() )
            return RefMetaDataTracker.EMPTY_TRACKER;
        else {
            List<RODRecordList> bindings = new ArrayList<RODRecordList>(states.size());

            for ( ReferenceOrderedDataState state: states )
                // todo -- warning, I removed the reference to the name from states
                bindings.add( state.iterator.seekForward(loc) );

            return new RefMetaDataTracker(bindings);
        }
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