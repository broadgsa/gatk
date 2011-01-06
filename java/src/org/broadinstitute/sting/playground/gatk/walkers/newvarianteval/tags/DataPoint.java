package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.tags;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface DataPoint {
    String description() default ""; // the description, optional
}
