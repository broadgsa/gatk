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

import org.broadinstitute.sting.utils.StingException;

import java.lang.reflect.Field;
import java.util.Collection;

/**
 * Describes the source field which defines a command-line argument.
 * A parsed-object version of the command-line argument will be
 * injected into an object containing this field.
 *
 * @author mhanna
 * @version 0.1
 */
public class ArgumentSource {
    /**
     * Class to which the field belongs.
     */
    public final Class clazz;

    /**
     * Field into which to inject command-line arguments.
     */
    public final Field field;

    /**
     * Descriptor for the argument.  Contains name, validation info, etc.
     */
    public final Argument descriptor;

    /**
     * Create a new command-line argument target.
     * @param clazz Class containing the argument.
     * @param field Field containing the argument.  Field must be annotated with 'Argument'.
     */
    public ArgumentSource( Class clazz, Field field ) {
        this.clazz = clazz;
        this.field = field;
        this.descriptor = field.getAnnotation(Argument.class);
        if( descriptor == null )
            throw new StingException("Cannot build out a command-line argument without a descriptor.");
    }

    /**
     * True if this argument source equals other.
     * @param other Another object, possibly an argument source, to test for equality.  Any object can
     *              be tested, but only instances of ArgumentSource will result in equals returning true.
     * @return True if this argument source matches other.  False otherwise.
     */
    @Override
    public boolean equals( Object other ) {
        if( other == null )
            return false;
        if( !(other instanceof ArgumentSource) )
            return false;

        ArgumentSource otherArgumentSource = (ArgumentSource)other;
        return this.clazz.equals(otherArgumentSource.clazz) && this.field.equals(otherArgumentSource.field);
    }

    /**
     * Returns an appropriate hash code for this argument source.
     * @return A uniformly distributed hashcode representing this argument source.
     */
    @Override
    public int hashCode() {
        return clazz.hashCode() ^ field.hashCode();
    }

    /**
     * Retrieves the full name of the argument, specifiable with the '--' prefix.  The full name can be
     * either specified explicitly with the fullName annotation parameter or implied by the field name.
     * @return full name of the argument.  Never null.
     */
    public String getFullName() {
        return descriptor.fullName().trim().length() > 0 ? descriptor.fullName().trim() : field.getName().toLowerCase();
    }

    /**
     * Retrieves the short name of the argument, specifiable with the '-' prefix.  The short name can
     * be specified or not; if left unspecified, no short name will be present.
     * @return short name of the argument.  Null if no short name exists.
     */
    public String getShortName() {
        return descriptor.shortName().trim().length() > 0 ? descriptor.shortName().trim() : null;
    }

    /**
     * Documentation for this argument.  Mandatory field.
     * @return Documentation for this argument.
     */
    public String getDoc() {
        return descriptor.doc();    
    }

    /**
     * Returns whether this field is required.  Note that flag fields are always forced to 'not required'.
     * @return True if the field is mandatory and not a boolean flag.  False otherwise.
     */
    public boolean isRequired() {
        return descriptor.required() && !isFlag();    
    }

    /**
     * Specifies other arguments which cannot be used in conjunction with tihs argument.  Comma-separated list.
     * @return A comma-separated list of exclusive arguments, or null if none are present.
     */
    public String getExclusiveOf() {
        return descriptor.exclusiveOf().trim().length() > 0 ? descriptor.exclusiveOf().trim() : null;    
    }

    /**
     * A regular expression which can be used for validation.
     * @return a JVM regex-compatible regular expression, or null to permit any possible value.
     */
    public String getValidationRegex() {
        return descriptor.validation().trim().length() > 0 ? descriptor.validation().trim() : null;    
    }

    /**
     * Set the value of the field in the passed object to <code>value</code>
     * @param targetInstance Instance in which to find the field.
     * @param value Value to which to set the field.
     */
    public void setValue( Object targetInstance, Object value ) {
        try {
            field.setAccessible(true);
            field.set(targetInstance, value);
        }
        catch( IllegalAccessException ex ) {
            //logger.fatal("processArgs: cannot convert field " + field.toString());
            throw new StingException("processArgs: Failed conversion " + ex.getMessage(), ex);
        }

    }

    /**
     * Returns true if the argument is a flag (a 0-valued argument).
     * @return True if argument is a flag; false otherwise.
     */
    public boolean isFlag() {
        return (field.getType() == Boolean.class) || (field.getType() == Boolean.TYPE);
    }

    /**
     * Can this argument support multiple values, or just one?
     * @return True if the argument supports multiple values.
     */
    public boolean isMultiValued() {
        Class argumentType = field.getType();
        return Collection.class.isAssignableFrom(argumentType) || field.getType().isArray();
    }


}
