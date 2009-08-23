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
     * Is this argument exclusive of other arguments?
     */
    public final String exclusiveOf;

    /**
     * Can we validate this regular expression?
     */
    public final String validation;

    /**
     * The target into which to inject arguments meeting this definition.
     */
    public final ArgumentSource source;

    /**
     * Creates a new argument definition.
     * @param source Source information for defining the argument.
     */
    public ArgumentDefinition( ArgumentSource source,
                               String fullName,
                               String shortName,
                               String doc,
                               boolean required,
                               String exclusiveOf,
                               String validation ) {
        this.source = source;

        this.fullName = fullName;
        this.shortName = shortName;
        this.doc = doc;
        this.required = required;
        this.exclusiveOf = exclusiveOf;
        this.validation = validation;
    }

    @Override
    public int hashCode() {
        int hashCode = fullName.hashCode();
        if(shortName != null) hashCode ^= shortName.hashCode();
        if(doc != null) hashCode ^= doc.hashCode();
        hashCode ^= Boolean.valueOf(required).hashCode();
        if(exclusiveOf != null) hashCode ^= exclusiveOf.hashCode();
        if(validation != null) hashCode ^= validation.hashCode();
        return hashCode;
    }

    public boolean equals( Object o ) {
        if( o == null )
            return false;
        if( !(o instanceof ArgumentDefinition) )
            return false;

        ArgumentDefinition other = (ArgumentDefinition)o;

        return Utils.equals(fullName,other.fullName) &&
               Utils.equals(shortName,other.shortName) &&
               Utils.equals(doc,other.doc) &&
               required == other.required &&
               Utils.equals(exclusiveOf,other.exclusiveOf) &&
               Utils.equals(validation,other.validation);
    }
}
