package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.tags;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface Analysis {
    String name() default ""; // its description, required
    String description(); // its description, required
}
