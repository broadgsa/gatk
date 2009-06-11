package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.apache.log4j.Logger;
/**
 * User: hanna
 * Date: May 12, 2009
 * Time: 11:24:42 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A queue of locus contexts.  Provides unidirectional seek.  Stripped down
 * implementation of java.util.Queue interface.
 */

public class CoveredLocusView extends LocusView {
    /**
     * Gets the position to which the last seek was requested.
     */
    private GenomeLoc seekPoint;

    /**
     * What's the context for the last locus accessed?
     * @param provider
     */
    private LocusContext nextLocusContext = null;

    private static Logger logger = Logger.getLogger(CoveredLocusView.class);

    /**
     * Create a new queue of locus contexts.
     * @param provider
     */
    public CoveredLocusView(ShardDataProvider provider) {
        super(provider);
    }

    public boolean hasNext() {
        return hasNextLocusContext();
    }

    public LocusContext next() {
        return nextLocusContext();
    }
}
