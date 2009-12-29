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

import java.lang.reflect.Field;
import java.util.Collection;
import java.util.List;

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
     * Type descriptor to use when parsing new argument types.
     */
    private final ArgumentTypeDescriptor typeDescriptor;

    /**
     * Create a new command-line argument target.
     * @param clazz Class containing the argument.
     * @param field Field containing the argument.  Field must be annotated with 'Argument'.
     */
    public ArgumentSource( Class clazz, Field field ) {
        this.clazz = clazz;
        this.field = field;
        this.typeDescriptor = ArgumentTypeDescriptor.create( field.getType() );
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
     * Generate a list of all argument definitions to which this argument source maps.
     * @return A non-null, non-empty list of argument definitions.
     */
    public List<ArgumentDefinition> createArgumentDefinitions() {
        return typeDescriptor.createArgumentDefinitions( this );
    }

    /**
     * This argument type descriptor wants to override any default value the user might have specified.
     * @return True if this descriptor wants to override any default the user specified.  False otherwise.
     */
    public boolean overridesDefault() {
        return typeDescriptor.overridesDefault();
    }

    /**
     * Provides the default value for the command-line argument.
     * @return Default value to load into the object.
     */
    public Object getDefault() {
        return typeDescriptor.getDefault();
    }

    /**
     * Parses the specified value based on the specified type.
     * @param source The type of value to be parsed.
     * @param values String representation of all values passed.
     * @return the parsed value of the object.
     */
    public Object parse( ArgumentSource source, ArgumentMatches values ) {
        Object value = null;
        if( !isFlag() )
            value = typeDescriptor.parse( source, values );
        else
            value = true;

        return value;
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

    /**
     * Gets a string representation of the argument source for debugging.
     * @return String representation of the argument source.
     */
    public String toString() {
        return clazz.getSimpleName() + ": " + field.getName();    
    }
}
