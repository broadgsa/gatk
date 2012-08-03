package org.broadinstitute.sting.gatk.walkers;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 8/2/12
 * Time: 1:58 PM
 * To change this template use File | Settings | File Templates.
 */

import java.lang.annotation.*;

/**
 * Indicates that program records should be removed from SAM headers by default for this walker
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface RemoveProgramRecords {
}
