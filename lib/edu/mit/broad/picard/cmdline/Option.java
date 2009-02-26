/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.cmdline;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Used to annotate which fields of a CommandLineProgram are options given at the command line.
 * If a command line call looks like "cmd option=foo x=y bar baz" the CommandLineProgram
 * would have annotations on fields to handle the values of option and x. All options
 * must be in the form name=value on the command line. The java type of the option
 * will be inferred from the type of the field or from the generic type of the collection
 * if this option is allowed more than once. The type must be an enum or
 * have a constructor with a single String parameter.
 *
 * @author Alec Wysoker
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
@Documented
public @interface Option {
	/** The name of the option as it would appear on the command line. */
    String shortName() default "";
    
    /** Text that appears for this option in text describing usage of the command line program. */
    String doc() default "";
    
    /**
     * If set to false, an exception will be thrown if the option is not specified.
     * If 2 options are mutually exclusive and both have optional=false it will be
     * interpreted as one or the other is required and an exception will only be thrown if
     * neither are specified. 
     */
    boolean optional() default false;
    
    /** 
     * Array of option names that cannot be used in conjunction with this one.
     * If 2 options are mutually exclusive and both have optional=false it will be
     * interpreted as one OR the other is required and an exception will only be thrown if
     * neither are specified. 
     */ 
    String[] mutex() default {};
    
    /** The minimum number of times that this option is required. */
    int minElements() default 0;
    
    /** The maximum number of times this option is allowed. */
    int maxElements() default Integer.MAX_VALUE;
}
