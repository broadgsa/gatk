package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.DownsampleType;

import java.lang.annotation.*;

/**
 * Specifies a method for downsampling the reads passed to a given
 * walker based on the input from that walker.
 *
 * @author hanna
 * @version 0.1
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface Downsample {
    DownsampleType by();
    int toCoverage() default -1;
    double toFraction() default -1.0F;
}