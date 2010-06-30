/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.commandline;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

import java.lang.annotation.Annotation;
import java.util.List;

/**
 * A specific argument definition.  Maps one-to-one with a field in some class.
 */
public class ArgumentDefinition {
    /**
     * Whether an argument is an input or an output.
     */
    public final ArgumentIOType ioType;

    /**
     * Full name of the argument.  Must have a value.
     */
    public final String fullName;

    /**
     * Short name of the argument.  Can be null.
     */
    public final String shortName;

    /**
     * Doc string for the argument.  Displayed in help.
     */
    public final String doc;

    /**
     * Is this argument required?
     */
    public final boolean required;

    /**
     * Is this argument a flag?  Users can't specify a value for a flag.
     */
    public final boolean isFlag;

    /**
     * Does this argument support multiple values (repeated "-arg value1 -arg value2"-style structures).
     */
    public final boolean isMultiValued;

    /**
     * Is this argument exclusive of other arguments?
     */
    public final String exclusiveOf;

    /**
     * Can we validate this regular expression?
     */
    public final String validation;

    /**
     * A list of valid options for this argument, if there is a compelling subset.
     */
    public final List<String> validOptions;

    /**
     * Creates a new argument definition.
     * @param ioType Whether the argument is an input or an output.
     * @param fullName Full name for this argument definition.
     * @param shortName Short name for this argument definition.
     * @param doc Doc string for this argument.
     * @param required Whether or not this argument is required.
     * @param isFlag Whether or not this argument should be treated as a flag.
     * @param isMultiValued Whether or not this argument supports multiple values.
     * @param exclusiveOf Whether this command line argument is mutually exclusive of other arguments.
     * @param validation A regular expression for command-line argument validation.
     * @param validOptions is there a particular list of options that's valid for this argument definition?  List them if so, otherwise set this to null. 
     */
    public ArgumentDefinition( ArgumentIOType ioType,
                               String fullName,
                               String shortName,
                               String doc,
                               boolean required,
                               boolean isFlag,
                               boolean isMultiValued,
                               String exclusiveOf,
                               String validation,
                               List<String> validOptions) {
        this.ioType = ioType;
        this.fullName = fullName;
        this.shortName = shortName;
        this.doc = doc;
        this.required = required;
        this.isFlag = isFlag;
        this.isMultiValued = isMultiValued;
        this.exclusiveOf = exclusiveOf;
        this.validation = validation;
        this.validOptions = validOptions;
    }

    /**
     * Creates a new argument definition.
     * @param annotation The annotation on the field.
     * @param defaultFullName Default full name for this argument definition.
     * @param defaultShortName Default short name for this argument definition.
     * @param isFlag Whether or not this argument should be treated as a flag.
     * @param isMultiValued Whether or not this argument supports multiple values.
     * @param validOptions is there a particular list of options that's valid for this argument definition?  List them if so, otherwise set this to null.
     */
    public ArgumentDefinition( Annotation annotation,
                               String defaultFullName,
                               String defaultShortName,
                               boolean isFlag,
                               boolean isMultiValued,
                               List<String> validOptions) {

        String fullName = (String)getValue(annotation, "fullName");
        String shortName = (String)getValue(annotation, "shortName");
        boolean isFullNameProvided = fullName.trim().length() > 0;
        boolean isShortNameProvided = shortName.trim().length() > 0;

        fullName = isFullNameProvided ? fullName.trim() : defaultFullName;

        // If the short name is provided, use that.  If the user hasn't provided any names at all, use
        // the default.  If somewhere in the middle, leave the short name blank.
        if( isShortNameProvided )
            shortName = shortName.trim();
        else if( !isFullNameProvided )
            shortName = defaultShortName;
        else
            shortName = null;

        this.ioType = getIOType(annotation);
        this.fullName = fullName;
        this.shortName = shortName;
        this.doc = getDoc(annotation);
        this.required = isRequired(annotation, isFlag);
        this.isFlag = isFlag;
        this.isMultiValued = isMultiValued;
        this.exclusiveOf = getExclusiveOf(annotation);
        this.validation = getValidationRegex(annotation);
        this.validOptions = validOptions;
    }
    
    /**
     * Creates a new argument definition.
     * @param annotation The annotation on the field.
     * @param fieldName Default full name for this argument definition.
     * @param isFlag Whether or not this argument should be treated as a flag.
     * @param isMultiValued Whether or not this argument supports multiple values.
     * @param validOptions is there a particular list of options that's valid for this argument definition?  List them if so, otherwise set this to null.
     */
    public ArgumentDefinition( Annotation annotation,
                               String fieldName,
                               boolean isFlag,
                               boolean isMultiValued,
                               List<String> validOptions) {
        this.ioType = getIOType(annotation);
        this.fullName = getFullName(annotation, fieldName);
        this.shortName = getShortName(annotation);
        this.doc = getDoc(annotation);
        this.required = isRequired(annotation, isFlag);
        this.isFlag = isFlag;
        this.isMultiValued = isMultiValued;
        this.exclusiveOf = getExclusiveOf(annotation);
        this.validation = getValidationRegex(annotation);
        this.validOptions = validOptions;
    }

    @Override
    public int hashCode() {
        int hashCode = fullName.hashCode();
        if(shortName != null) hashCode ^= shortName.hashCode();
        return hashCode;
    }

    public boolean equals( Object o ) {
        if( o == null )
            return false;
        if( !(o instanceof ArgumentDefinition) )
            return false;

        ArgumentDefinition other = (ArgumentDefinition)o;

        return Utils.equals(fullName,other.fullName) &&
               Utils.equals(shortName,other.shortName);
    }

    /**
     * Returns the ArgumentIOType for the annotation.
     * @param annotation @Input or @Output
     * @return ArgumentIOType.Input, Output, or Unknown
     */
    public static ArgumentIOType getIOType(Annotation annotation) {
        if (annotation instanceof Input) return ArgumentIOType.INPUT;
        if (annotation instanceof Output) return ArgumentIOType.OUTPUT;
        return ArgumentIOType.UNKNOWN;
    }
    
    /**
     * A hack to get around the fact that Java doesn't like inheritance in Annotations.
     * @param annotation to run the method on
     * @param method the method to invoke
     * @return the return value of the method
     */
    private static Object getValue(Annotation annotation, String method) {
        try {
            return annotation.getClass().getMethod(method).invoke(annotation);
        } catch (Exception e) {
            throw new StingException("Unable to access method " + method + " on annotation " + annotation.getClass(), e);
        }
    }

    /**
     * Retrieves the full name of the argument, specifiable with the '--' prefix.  The full name can be
     * either specified explicitly with the fullName annotation parameter or implied by the field name.
     * @param annotation Original field annotation.
     * @param fieldName Original field name.
     * @return full name of the argument.  Never null.
     */
    private static String getFullName( Annotation annotation, String fieldName ) {
        String fullName = (String)getValue(annotation, "fullName");
        return fullName.trim().length() > 0 ? fullName.trim() : fieldName.toLowerCase();
    }

    /**
     * Retrieves the short name of the argument, specifiable with the '-' prefix.  The short name can
     * be specified or not; if left unspecified, no short name will be present.
     * @param annotation Original field annotation.
     * @return short name of the argument.  Null if no short name exists.
     */
    private static String getShortName( Annotation annotation ) {
        String shortName = (String)getValue(annotation, "shortName");
        return shortName.trim().length() > 0 ? shortName.trim() : null;
    }

    /**
     * Documentation for this argument.  Mandatory field.
     * @param annotation Original field annotation.
     * @return Documentation for this argument.
     */
    private static String getDoc( Annotation annotation ) {
        return (String)getValue(annotation, "doc");
    }

    /**
     * Returns whether this field is required.  Note that flag fields are always forced to 'not required'.
     * @param annotation Original field annotation.
     * @param isFlag True if the field is a flag.
     * @return True if the field is mandatory and not a boolean flag.  False otherwise.
     */
    private static boolean isRequired( Annotation annotation, boolean isFlag ) {
        boolean required = (Boolean)getValue(annotation, "required");
        return required && !isFlag;
    }

    /**
     * Specifies other arguments which cannot be used in conjunction with this argument.  Comma-separated list.
     * @param annotation Original field annotation.
     * @return A comma-separated list of exclusive arguments, or null if none are present.
     */
    private static String getExclusiveOf( Annotation annotation ) {
        String exclusiveOf = (String)getValue(annotation, "exclusiveOf");
        return exclusiveOf.trim().length() > 0 ? exclusiveOf.trim() : null;
    }

    /**
     * A regular expression which can be used for validation.
     * @param annotation Original field annotation.
     * @return a JVM regex-compatible regular expression, or null to permit any possible value.
     */
    private static String getValidationRegex( Annotation annotation ) {
        String validation = (String)getValue(annotation, "validation");
        return validation.trim().length() > 0 ? validation.trim() : null;
    }
}
