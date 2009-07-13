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
public abstract class FieldParser {

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(FieldParser.class);

    /**
     * Name of the field which should be parsed by this argument.
     */
    protected final String fieldName;

    public static FieldParser create( Field field ) {
        Class type = field.getType();

        if( Collection.class.isAssignableFrom(type) || type.isArray() )
            return new JRECompoundFieldParser( field );
        else
            return new JRESimpleFieldParser( field );
    }

    static FieldParser create( String fieldName, Class type ) {
        if( Collection.class.isAssignableFrom(type) || type.isArray() )
            return new JRECompoundFieldParser( fieldName, type );
        else
            return new JRESimpleFieldParser( fieldName, type );        
    }

    protected FieldParser( Field field ) {
        fieldName = field.toString();
    }

    protected FieldParser( String fieldName ) {
        this.fieldName = fieldName;
    }
    
    public abstract Object parse( String... values );
}

class JRESimpleFieldParser extends FieldParser {
    private final Class type;

    public JRESimpleFieldParser( Field field ) {
        super( field );
        this.type = field.getType();
    }

    public JRESimpleFieldParser( String fieldName, Class type ) {
        super( fieldName );
        this.type = type;
    }

    public Object parse( String... values ) {
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
            logger.fatal("ArgumentParser: NoSuchMethodException: cannot convert field " + fieldName);
            throw new StingException("constructFromString:NoSuchMethodException: Failed conversion " + e.getMessage());
        } catch (IllegalAccessException e) {
            logger.fatal("ArgumentParser: IllegalAccessException: cannot convert field " + fieldName);
            throw new StingException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
        } catch (InvocationTargetException e) {
            logger.fatal("ArgumentParser: InvocationTargetException: cannot convert field " + fieldName);
            throw new StingException("constructFromString:InvocationTargetException: Failed conversion " + e.getMessage());
        } catch (InstantiationException e) {
            logger.fatal("ArgumentParser: InstantiationException: cannot convert field " + fieldName);
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

class JRECompoundFieldParser extends FieldParser {
    private final Class type;
    private final Class componentType;
    private final FieldParser componentArgumentParser;

    public JRECompoundFieldParser( Field field ) {
        super( field );
        Class candidateType = field.getType();

        if( Collection.class.isAssignableFrom(candidateType) ) {

            // If this is a generic interface, pick a concrete implementation to create and pass back.
            // Because of type erasure, don't worry about creating one of exactly the correct type.
            if( Modifier.isInterface(candidateType.getModifiers()) || Modifier.isAbstract(candidateType.getModifiers()) )
            {
                if( java.util.List.class.isAssignableFrom(candidateType) ) candidateType = ArrayList.class;
                else if( java.util.Queue.class.isAssignableFrom(candidateType) ) candidateType = java.util.ArrayDeque.class;
                else if( java.util.Set.class.isAssignableFrom(candidateType) ) candidateType = java.util.TreeSet.class;
            }

            this.type = candidateType;

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
        else if( candidateType.isArray() ) {
            this.type = candidateType;
            this.componentType = candidateType.getComponentType();
        }
        else
            throw new StingException("Unsupported compound argument type: " + candidateType);

        componentArgumentParser = FieldParser.create( fieldName, componentType );
    }

    public JRECompoundFieldParser( String fieldName, Class type ) {
        super(fieldName);

        this.type = type;

        if( Collection.class.isAssignableFrom(type) ) {
            this.componentType = String.class;
        }
        else if( type.isArray() ) {
            this.componentType = type.getComponentType();
        }
        else
            throw new StingException("Unsupported compound argument type: " + type);

        componentArgumentParser = FieldParser.create( fieldName, componentType );        
    }

    @Override
    public Object parse( String... values ) 
    {
        if( Collection.class.isAssignableFrom(type) ) {
            Collection collection = null;
            try {
                collection = (Collection)type.newInstance();
            }
            catch (InstantiationException e) {
                logger.fatal("ArgumentParser: InstantiationException: cannot convert field " + fieldName);
                throw new StingException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
            }
            catch (IllegalAccessException e) {
                logger.fatal("ArgumentParser: IllegalAccessException: cannot convert field " + fieldName);
                throw new StingException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
            }

            for( String value: values )
                collection.add( componentArgumentParser.parse(value) );

            return collection;
        }
        else if( type.isArray() ) {
            Object arr = Array.newInstance(componentType,values.length);

            for( int i = 0; i < values.length; i++ )
                Array.set( arr,i,componentArgumentParser.parse(values[i]));
            
            return arr;
        }
        else
            throw new StingException("Unsupported compound argument type: " + type);
    }
}

