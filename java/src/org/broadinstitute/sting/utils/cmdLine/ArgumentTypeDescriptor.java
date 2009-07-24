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
import org.apache.log4j.Logger;

import java.lang.reflect.*;
import java.util.*;

/**
 * An factory capable of providing parsers that can parse any type
 * of supported command-line argument.
 *
 * @author mhanna
 * @version 0.1
 */
public abstract class ArgumentTypeDescriptor {

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(ArgumentTypeDescriptor.class);
    
    /**
     * Class reference to the different types of descriptors that the create method can create.
     */
    private static List<ArgumentTypeDescriptor> descriptors = Arrays.asList( new SimpleArgumentTypeDescriptor(),                                                                                
                                                                             new CompoundArgumentTypeDescriptor() );

    public static ArgumentTypeDescriptor create( Class type ) {
        for( ArgumentTypeDescriptor descriptor: descriptors ) {
            if( descriptor.supports(type) )
                return descriptor;
        }
        throw new StingException("Can't process command-line arguments of type: " + type.getName());
    }

    public abstract boolean supports( Class type );
    
    public abstract Object parse( Field field, Class type, String... values );
}

class SimpleArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    @Override
    public boolean supports( Class type ) {
        if( type.isPrimitive() ) return true;
        if( type.isEnum() ) return true;
        if( primitiveToWrapperMap.containsValue(type) ) return true;

        try {
            type.getConstructor(String.class);
            return true;
        }
        catch( Exception ex ) {
            // An exception thrown above means that the String constructor either doesn't
            // exist or can't be accessed.  In either case, this descriptor doesn't support this type.
            return false;
        }
    }

    @Override
    public Object parse( Field field, Class type, String... values ) {
        if( values.length > 1 )
            throw new StingException("Simple argument parser is unable to parse multiple arguments.");

        String value = values[0];

        // lets go through the types we support
        try {
            if (type.isPrimitive()) {
                Method valueOf = primitiveToWrapperMap.get(type).getMethod("valueOf",String.class);
                return valueOf.invoke(null,value.trim());
            } else if (type.isEnum()) {
                return Enum.valueOf(type,value.toUpperCase().trim());
            } else {
                Constructor ctor = type.getConstructor(String.class);
                return ctor.newInstance(value);
            }
        }
        catch (NoSuchMethodException e) {
            throw new StingException("constructFromString:NoSuchMethodException: Failed conversion " + e.getMessage());
        } catch (IllegalAccessException e) {
            throw new StingException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
        } catch (InvocationTargetException e) {
            throw new StingException("constructFromString:InvocationTargetException: Failed conversion " + e.getMessage());
        } catch (InstantiationException e) {
            throw new StingException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
        }

    }

    /**
     * A mapping of the primitive types to their associated wrapper classes.  Is there really no way to infer
     * this association available in the JRE?
     */
    private static Map<Class,Class> primitiveToWrapperMap = new HashMap<Class,Class>() {
        {
            put( Boolean.TYPE, Boolean.class );
            put( Character.TYPE, Character.class );
            put( Byte.TYPE, Byte.class );
            put( Short.TYPE, Short.class );
            put( Integer.TYPE, Integer.class );
            put( Long.TYPE, Long.class );
            put( Float.TYPE, Float.class );
            put( Double.TYPE, Double.class );
        }
    };    
}

class CompoundArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    @Override
    public boolean supports( Class type ) {
        return ( Collection.class.isAssignableFrom(type) || type.isArray() );
    }
    
    @Override
    public Object parse( Field field, Class type, String... values )
    {
        Class componentType = null;
        ArgumentTypeDescriptor componentArgumentParser;        

        if( Collection.class.isAssignableFrom(type) ) {

            // If this is a generic interface, pick a concrete implementation to create and pass back.
            // Because of type erasure, don't worry about creating one of exactly the correct type.
            if( Modifier.isInterface(type.getModifiers()) || Modifier.isAbstract(type.getModifiers()) )
            {
                if( java.util.List.class.isAssignableFrom(type) ) type = ArrayList.class;
                else if( java.util.Queue.class.isAssignableFrom(type) ) type = java.util.ArrayDeque.class;
                else if( java.util.Set.class.isAssignableFrom(type) ) type = java.util.TreeSet.class;
            }

            // If this is a parameterized collection, find the contained type.  If blow up if only one type exists.
            if( field.getGenericType() instanceof ParameterizedType) {
                ParameterizedType parameterizedType = (ParameterizedType)field.getGenericType();
                if( parameterizedType.getActualTypeArguments().length > 1 )
                    throw new IllegalArgumentException("Unable to determine collection type of field: " + field.toString());
                componentType = (Class)parameterizedType.getActualTypeArguments()[0];
            }
            else
                componentType = String.class;
        }
        else if( type.isArray() ) {
            componentType = type.getComponentType();
        }
        else
            throw new StingException("Unsupported compound argument type: " + type);

        componentArgumentParser = ArgumentTypeDescriptor.create( componentType );

        if( Collection.class.isAssignableFrom(type) ) {
            Collection collection = null;
            try {
                collection = (Collection)type.newInstance();
            }
            catch (InstantiationException e) {
                logger.fatal("ArgumentParser: InstantiationException: cannot convert field " + field.getName());
                throw new StingException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
            }
            catch (IllegalAccessException e) {
                logger.fatal("ArgumentParser: IllegalAccessException: cannot convert field " + field.getName());
                throw new StingException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
            }

            for( String value: values )
                collection.add( componentArgumentParser.parse(field,componentType,value) );

            return collection;
        }
        else if( type.isArray() ) {
            Object arr = Array.newInstance(componentType,values.length);

            for( int i = 0; i < values.length; i++ )
                Array.set( arr,i,componentArgumentParser.parse(field,componentType,values[i]));
            
            return arr;
        }
        else
            throw new StingException("Unsupported compound argument type: " + type);
    }
}

