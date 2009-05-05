package org.broadinstitute.sting.utils.cmdLine;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 24, 2009
 * Time: 11:11:36 AM
 */
/**
 * Annotates fields in objects that should be used as command-line arguments.
 * Any field annotated with @Argument can appear as a command-line parameter. 
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.FIELD})
public @interface Argument {
    /**
     * The full name of the command-line argument.  Full names should be
     * prefixed on the command-line with a double dash (--).
     * @return Selected full name, or "" to use the default.
     */
    String fullName() default "";

    /**
     * Specified short name of the command.  Short names should be prefixed
     * with a single dash.  Argument values can directly abut single-char
     * short names or be separated from them by a space.
     * @return Selected short name, or "" for none.
     */
    String shortName() default "";

    /**
     * Documentation for the command-line argument.  Should appear when the
     * --help argument is specified. 
     * @return Doc string associated with this command-line argument.
     */
    String doc();

    /**
     * Is this command-line argument required.  The application should exit
     * printing help if this command-line argument is not specified.
     * @return True if the argument is required.  False otherwise.
     */
    boolean required() default true;

    /**
     * Should this command-line argument be exclusive of others.  Should be
     * a comma-separated list of names of arguments of which this should be
     * independent. 
     * @return A comma-separated string listing other arguments of which this
     *         argument should be independent.
     */
    String exclusive() default "";

    /**
     * What is the default value for this argument type.
     * @return Default value of this argument type.
     */
    String defaultValue() default "";
}
