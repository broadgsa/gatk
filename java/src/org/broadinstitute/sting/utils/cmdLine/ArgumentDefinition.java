/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.Utils;

import java.util.List;
import java.util.Collections;

/**
 * A specific argument definition.  Maps one-to-one with a field in some class.
 */
public class ArgumentDefinition {
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
    public ArgumentDefinition( String fullName,
                               String shortName,
                               String doc,
                               boolean required,
                               boolean isFlag,
                               boolean isMultiValued,
                               String exclusiveOf,
                               String validation,
                               List<String> validOptions) {
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
}
