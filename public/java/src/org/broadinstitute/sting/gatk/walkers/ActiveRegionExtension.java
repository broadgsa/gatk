package org.broadinstitute.sting.gatk.walkers;

import java.lang.annotation.Documented;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Describes the size of the buffer region that is added to each active region when pulling in covered reads.
 * User: rpoplin
 * Date: 1/18/12
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)

public @interface ActiveRegionExtension {
    public int extension() default 0;
    public int maxRegion() default 1500;
}
