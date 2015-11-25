/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.commandline;

import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

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
     * The class of the argument.
     */
    public final Class argumentType;

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
     * Must this argument be specified on the command-line?  Note that there's a
     * critical difference between the meaning of a required argument from the
     * perspective of the argument source and the perspective of the argument
     * definition: the argument source's required field indicates that the field
     * should somehow be populated by the GATK (and fail if there's an error).
     * The ArgumentDefinition required element means that the required element
     * must be specified on the command-line.
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
     * The class of the componentType.  Not used for scalars.
     */
    public final Class componentType;

    /**
     * Is this argument hidden from the help system?
     */
    public final boolean isHidden;

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
     * @param argumentType The class of the field.
     * @param fullName Full name for this argument definition.
     * @param shortName Short name for this argument definition.
     * @param doc Doc string for this argument.
     * @param required Whether or not this argument is required.
     * @param isFlag Whether or not this argument should be treated as a flag.
     * @param isMultiValued Whether or not this argument supports multiple values.
     * @param isHidden Whether or not this argument should be hidden from the command-line argument system.
     * @param componentType For multivalued arguments the type of the components.
     * @param exclusiveOf Whether this command line argument is mutually exclusive of other arguments.
     * @param validation A regular expression for command-line argument validation.
     * @param validOptions is there a particular list of options that's valid for this argument definition?  List them if so, otherwise set this to null. 
     */
    public ArgumentDefinition( ArgumentIOType ioType,
                               Class argumentType,
                               String fullName,
                               String shortName,
                               String doc,
                               boolean required,
                               boolean isFlag,
                               boolean isMultiValued,
                               boolean isHidden,
                               Class componentType,
                               String exclusiveOf,
                               String validation,
                               List<String> validOptions) {
        this.ioType = ioType;
        this.argumentType = argumentType;
        this.fullName = fullName;
        this.shortName = shortName;
        this.doc = doc;
        this.required = required;
        this.isFlag = isFlag;
        this.isMultiValued = isMultiValued;
        this.isHidden = isHidden;
        this.componentType = componentType;
        this.exclusiveOf = exclusiveOf;
        this.validation = validation;
        this.validOptions = validOptions;

        validateName(shortName);
        validateName(fullName);
    }

    /**
     * Creates a new argument definition.
     * @param annotation The annotation on the field.
     * @param argumentType The class of the field.
     * @param defaultFullName Default full name for this argument definition.
     * @param defaultShortName Default short name for this argument definition.
     * @param isFlag Whether or not this argument should be treated as a flag.
     * @param isMultiValued Whether or not this argument supports multiple values.
     * @param componentType For multivalued arguments the type of the components.
     * @param isHidden Whether or not this argument should be hidden from the command-line argument system.
     * @param validOptions is there a particular list of options that's valid for this argument definition?  List them if so, otherwise set this to null.
     */
    public ArgumentDefinition( Annotation annotation,
                               ArgumentIOType ioType,
                               Class argumentType,
                               String defaultFullName,
                               String defaultShortName,
                               String doc,
                               boolean isRequired,
                               boolean isFlag,
                               boolean isMultiValued,
                               boolean isHidden,
                               Class componentType,
                               String exclusiveOf,
                               String validation,
                               List<String> validOptions) {

        String fullName = (String)CommandLineUtils.getValue(annotation, "fullName");
        String shortName = (String)CommandLineUtils.getValue(annotation, "shortName");
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

        validateName(shortName);
        validateName(fullName);

        this.ioType = ioType;
        this.argumentType = argumentType;
        this.fullName = fullName;
        this.shortName = shortName;
        this.doc = doc;
        this.required = isRequired;
        this.isFlag = isFlag;
        this.isMultiValued = isMultiValued;
        this.isHidden = isHidden;
        this.componentType = componentType;
        this.exclusiveOf = exclusiveOf;
        this.validation = validation;
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
     * Retrieves the full name of the argument, specifiable with the '--' prefix.  The full name can be
     * either specified explicitly with the fullName annotation parameter or implied by the field name.
     * @param annotation Original field annotation.
     * @param fieldName Original field name.
     * @return full name of the argument.  Never null.
     */
    public static String getFullName( Annotation annotation, String fieldName ) {
        String fullName = (String)CommandLineUtils.getValue(annotation, "fullName");
        return fullName.trim().length() > 0 ? fullName.trim() : fieldName.toLowerCase();
    }

    /**
     * Retrieves the short name of the argument, specifiable with the '-' prefix.  The short name can
     * be specified or not; if left unspecified, no short name will be present.
     * @param annotation Original field annotation.
     * @return short name of the argument.  Null if no short name exists.
     */
    public static String getShortName( Annotation annotation ) {
        String shortName = (String)CommandLineUtils.getValue(annotation, "shortName");
        return shortName.trim().length() > 0 ? shortName.trim() : null;
    }

    /**
     * Documentation for this argument.  Mandatory field.
     * @param annotation Original field annotation.
     * @return Documentation for this argument.
     */
    public static String getDoc( Annotation annotation ) {
        return (String)CommandLineUtils.getValue(annotation, "doc");
    }

    /**
     * Specifies other arguments which cannot be used in conjunction with this argument.  Comma-separated list.
     * @param annotation Original field annotation.
     * @return A comma-separated list of exclusive arguments, or null if none are present.
     */
    public static String getExclusiveOf( Annotation annotation ) {
        String exclusiveOf = (String)CommandLineUtils.getValue(annotation, "exclusiveOf");
        return exclusiveOf.trim().length() > 0 ? exclusiveOf.trim() : null;
    }

    /**
     * A regular expression which can be used for validation.
     * @param annotation Original field annotation.
     * @return a JVM regex-compatible regular expression, or null to permit any possible value.
     */
    public static String getValidationRegex( Annotation annotation ) {
        String validation = (String)CommandLineUtils.getValue(annotation, "validation");
        return validation.trim().length() > 0 ? validation.trim() : null;
    }

    /**
     * Make sure the argument's name is valid
     *
     * @param name
     */
    private void validateName(final String name) {
        if ( name != null && name.startsWith("-") )
            throw new ReviewedGATKException("Invalid argument definition: " + name + " begins with a -");
    }
}
