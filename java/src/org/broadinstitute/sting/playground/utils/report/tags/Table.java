package org.broadinstitute.sting.playground.utils.report.tags;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * @author aaron
 *         <p/>
 *         Annotation Table
 *         <p/>
 *         tells the report system to make a table out of the annotated field.  The field
 *         can be of type List<>, Map<>, Vector<>, Collection<> or primative array.
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface Table {
    String name();              // the name
    String header() default ""; // any text representing the table header
    String description();       // a description of the table
    int columns() default -1;   // the number of columns to divide the data into
}
