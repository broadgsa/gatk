package org.broadinstitute.sting.gatk.walkers;

import java.lang.annotation.*;
/**
 * User: hanna
 * Date: May 19, 2009
 * Time: 10:05:01 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Determines what data sources are allowed by a given walker.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface Allows {
    DataSource[] value();
}
