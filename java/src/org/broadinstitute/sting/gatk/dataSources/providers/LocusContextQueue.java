package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;
/**
 * User: hanna
 * Date: May 13, 2009
 * Time: 3:30:16 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A queue of locus context entries.
 */

public interface LocusContextQueue {
    /**
     * Get the locus context at the given position.
     * @return Locus context, or null if no locus context exists at this position.
     */
    LocusContext peek();

    /**
     * Seek to the given point the queue of locus contexts.
     * @param target Target base pair to which to seek.  Must be a single base pair.
     * @return an instance of itself for parameter chaining.
     */
    public LocusContextQueue seek(GenomeLoc target);
}
