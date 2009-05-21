package org.broadinstitute.sting.gatk.dataSources.providers;
/**
 * User: hanna
 * Date: May 21, 2009
 * Time: 3:14:56 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Represents a view into given data.
 */
public interface View {
    /**
     * Inform this view that the data provided to it no longer exists.
     */
    public void close();
}
