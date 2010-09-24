package org.broadinstitute.sting.utils.report.tags;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * @author aaron
 *         <p/>
 *         Annotation Param
 *         <p/>
 *         a description annotation for a parameter; a variable used as input to the
 *         analysis, but not (nessasarally) an output.  Some formats will store this
 *         information in comments, others will not include it.
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface Param {
    String name() default "";   // the name, defaulted to the variable name
    String description();       // the description of the parameter
}
