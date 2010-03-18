package org.broadinstitute.sting.playground.utils.report.tags;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * @author aaron
 *         <p/>
 *         Annotation Analysis
 *         <p/>
 *         the main annotation for analysis objects in the report system
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface Analysis {
    String name() default "";               // the name of the analysis
    String description();                   // its description, required
    String version() default "";            // the version, not always used
}
