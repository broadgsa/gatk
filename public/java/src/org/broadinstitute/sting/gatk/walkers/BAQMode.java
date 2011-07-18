package org.broadinstitute.sting.gatk.walkers;

import java.lang.annotation.*;

/**
 * User: hanna
 * Date: May 14, 2009
 * Time: 1:51:22 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Allows the walker to indicate what type of data it wants to consume.
 */

@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface BAQMode {
    public abstract org.broadinstitute.sting.utils.baq.BAQ.QualityMode QualityMode() default org.broadinstitute.sting.utils.baq.BAQ.QualityMode.OVERWRITE_QUALS;
    public abstract org.broadinstitute.sting.utils.baq.BAQ.ApplicationTime ApplicationTime() default org.broadinstitute.sting.utils.baq.BAQ.ApplicationTime.ON_INPUT;
}