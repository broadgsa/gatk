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

package org.broadinstitute.gatk.engine.traversals;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.ReadMetrics;
import org.broadinstitute.gatk.engine.datasources.providers.ShardDataProvider;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.progressmeter.ProgressMeter;

public abstract class TraversalEngine<M,T,WalkerType extends Walker<M,T>,ProviderType extends ShardDataProvider> {
    /** our log, which we want to capture anything from this class */
    protected static final Logger logger = Logger.getLogger(TraversalEngine.class);

    protected GenomeAnalysisEngine engine;
    private ProgressMeter progressMeter;

    // ----------------------------------------------------------------------------------------------------
    //
    // ABSTRACT METHODS
    //
    // ----------------------------------------------------------------------------------------------------

    /**
     * Gets the named traversal type associated with the given traversal, such as loci, reads, etc.
     *
     * @return A user-friendly name for the given traversal type.
     */
    public abstract String getTraversalUnits();

    /**
     * this method must be implemented by all traversal engines
     *
     * @param walker       the walker to run with
     * @param dataProvider the data provider that generates data given the shard
     * @param sum          the accumulator
     *
     * @return an object of the reduce type
     */
    public abstract T traverse(WalkerType walker,
                               ProviderType dataProvider,
                               T sum);

    /**
     * Initialize the traversal engine.  After this point traversals can be run over the data
     *
     * @param engine GenomeAnalysisEngine for this traversal
     * @param progressMeter An optional (null == optional) meter to track our progress
     */
    public void initialize(final GenomeAnalysisEngine engine, final Walker walker, final ProgressMeter progressMeter) {
        if ( engine == null )
            throw new ReviewedGATKException("BUG: GenomeAnalysisEngine cannot be null!");

        this.engine = engine;
        this.progressMeter = progressMeter;
    }

    /**
     * For testing only.  Does not initialize the progress meter
     *
     * @param engine
     */
    protected void initialize(final GenomeAnalysisEngine engine, final Walker walker) {
        initialize(engine, walker, null);
    }

    /**
     * Called by the MicroScheduler when all work is done and the GATK is shutting down.
     *
     * To be used by subclasses that need to free up resources (such as threads)
     */
    public void shutdown() {
        // by default there's nothing to do
    }

    /**
     * Update the cumulative traversal metrics according to the data in this shard
     *
     * @param singleTraverseMetrics read metrics object containing the information about a single shard's worth
     *                              of data processing
     */
    public void updateCumulativeMetrics(final ReadMetrics singleTraverseMetrics) {
        engine.getCumulativeMetrics().incrementMetrics(singleTraverseMetrics);
    }

    /**
     * Forward request to notifyOfProgress
     *
     * Assumes that one cycle has been completed
     *
     * @param loc  the location
     */
    public void printProgress(final GenomeLoc loc) {
        if ( progressMeter != null )
            progressMeter.notifyOfProgress(loc, engine.getCumulativeMetrics().getNumIterations());
    }
}

