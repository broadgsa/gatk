package org.broadinstitute.sting.gatk.walkers;
/**
 * User: hanna
 * Date: May 14, 2009
 * Time: 2:12:33 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Allow user to choose between a number of different data sources.
 */
public enum DataSource {
    READS,
    REFERENCE,
    REFERENCE_BASES,         // Do I actually need the reference bases passed to the walker?
    REFERENCE_ORDERED_DATA
}
