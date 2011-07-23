package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface DataPoint {
    String description() default ""; // the description, optional
}
