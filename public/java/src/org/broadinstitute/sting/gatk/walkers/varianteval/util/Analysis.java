package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface Analysis {
    String name() default ""; // its description, required
    String description(); // its description, required
}
