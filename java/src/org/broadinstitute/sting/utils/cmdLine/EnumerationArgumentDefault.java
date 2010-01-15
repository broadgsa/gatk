package org.broadinstitute.sting.utils.cmdLine;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * @author aaron
 * <p/>
 * Annotation EnumerationArgumentDefault
 * <p/>
 * Allows the default argument value to be set for an enum; this allows us to treat enums as
 * booleans on the command line. I.e.
 *
 * if we're using an enum Shape,
 *
 * enum shape {
 *  SQUARE,
 *  CIRCLE,
 *  @EnumerationArgumentDefault
 *  TRIANGLE
 * }
 *
 * and a command line option -shape, the EnumerationArgumentDefault would allow you to say:
 * -shape
 * or
 * -shape TRIANGLE
 *
 * would get -shape set to TRIANGLE, where:
 *
 * -shape SQUARE
 *
 * would set shape to SQUARE
 *
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
public @interface EnumerationArgumentDefault {
}
