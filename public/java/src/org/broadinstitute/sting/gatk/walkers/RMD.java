package org.broadinstitute.sting.gatk.walkers;

import org.broad.tribble.Feature;

import java.lang.annotation.Documented;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
/**
 * User: hanna
 * Date: May 19, 2009
 * Time: 1:34:15 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A data type representing reference-ordered data.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
public @interface RMD {
    String name();    
    Class type() default Feature.class;
}
