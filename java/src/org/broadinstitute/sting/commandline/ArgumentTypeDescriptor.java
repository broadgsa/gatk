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
import org.apache.log4j.Logger;

import java.lang.reflect.*;
import java.util.*;

/**
 * An descriptor capable of providing parsers that can parse any type
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
     * The type of set used must be ordered (but not necessarily sorted).
     */
    private static Set<ArgumentTypeDescriptor> descriptors = new LinkedHashSet<ArgumentTypeDescriptor>( Arrays.asList(new SimpleArgumentTypeDescriptor(),
                                                                                                                      new CompoundArgumentTypeDescriptor()) );

    /**
     * Adds new, user defined descriptors to the head of the descriptor list.
     * @param argumentTypeDescriptors New descriptors to add.  List can be empty, but should not be null.
     */
    public static void addDescriptors( Collection<ArgumentTypeDescriptor> argumentTypeDescriptors ) {
        // We care about ordering; newly added descriptors should have priority over stock descriptors.
        // Enforce this by creating a new *ordered* set, adding the new descriptors, then adding the old descriptors.
        Set<ArgumentTypeDescriptor> allDescriptors = new LinkedHashSet<ArgumentTypeDescriptor>();
        allDescriptors.addAll( argumentTypeDescriptors );
        allDescriptors.addAll( descriptors );
        descriptors = allDescriptors;
    }

    /**
     * Fetch the given descriptor from the descriptor repository.
     * @param type Class for which to specify a descriptor.
     * @return descriptor for the given type.
     */
    public static ArgumentTypeDescriptor create( Class type ) {
        for( ArgumentTypeDescriptor descriptor: descriptors ) {
            if( descriptor.supports(type) )
                return descriptor;
        }
        throw new StingException("Can't process command-line arguments of type: " + type.getName());
    }

    /**
     * Does this descriptor support classes of the given type?
     * @param type The type to check.
     * @return true if this descriptor supports the given type, false otherwise.
     */
    public abstract boolean supports( Class type );

    /**
     * This argument type descriptor wants to override any default value the user might have specified. 
     * @return True if this descriptor wants to override any default the user specified.  False otherwise.
     */
    public boolean overridesDefault() {
        return false;
    }

    /**
     * Provides the default value for the command-line argument.
     * @return Default value to load into the object.
     */
    public Object getDefault() {
        throw new UnsupportedOperationException(String.format("Type descriptor %s cannot override default value of command-line argument",this.getClass()));   
    }

    /**
     * Given the given argument source and attributes, synthesize argument definitions for command-line arguments.
     * @param source Source class and field for the given argument.
     * @return A list of command-line argument definitions supporting this field.
     */
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source ) {
        return Collections.singletonList(createDefaultArgumentDefinition(source));
    }

    public Object parse( ArgumentSource source, ArgumentMatches matches ) {
        return parse( source, source.field.getType(), matches );
    }

    /**
     * By default, argument sources create argument definitions with a set of default values.
     * Use this method to create the one simple argument definition.
     * @param source argument source for which to create a default definition.
     * @return The default definition for this argument source.
     */
    protected ArgumentDefinition createDefaultArgumentDefinition( ArgumentSource source ) {
        return new ArgumentDefinition( getFullName(source),
                                       getShortName(source),
                                       getDoc(source),
                                       isRequired(source),
                                       source.isFlag(),
                                       source.isMultiValued(),
                                       getExclusiveOf(source),
                                       getValidationRegex(source),
                                       getValidOptions(source) );
    }

    protected abstract Object parse( ArgumentSource source, Class type, ArgumentMatches matches );

    /**
     * Retrieves the full name of the argument, specifiable with the '--' prefix.  The full name can be
     * either specified explicitly with the fullName annotation parameter or implied by the field name.
     * @param source Original field specifying command-line arguments.
     * @return full name of the argument.  Never null.
     */
    protected String getFullName( ArgumentSource source ) {
        Argument description = getArgumentDescription(source);
        return description.fullName().trim().length() > 0 ? description.fullName().trim() : source.field.getName().toLowerCase();
    }

    /**
     * Retrieves the short name of the argument, specifiable with the '-' prefix.  The short name can
     * be specified or not; if left unspecified, no short name will be present.
     * @param source Original field specifying command-line arguments.
     * @return short name of the argument.  Null if no short name exists.
     */
    protected String getShortName( ArgumentSource source ) {
        Argument description = getArgumentDescription(source);
        return description.shortName().trim().length() > 0 ? description.shortName().trim() : null;
    }

    /**
     * Documentation for this argument.  Mandatory field.
     * @param source Original field specifying command-line arguments.
     * @return Documentation for this argument.
     */
    protected String getDoc( ArgumentSource source ) {
        Argument description = getArgumentDescription(source);
        return description.doc();
    }

    /**
     * Returns whether this field is required.  Note that flag fields are always forced to 'not required'.
     * @param source Original field specifying command-line arguments.
     * @return True if the field is mandatory and not a boolean flag.  False otherwise.
     */
    protected boolean isRequired( ArgumentSource source ) {
        Argument description = getArgumentDescription(source);
        return description.required() && !source.isFlag();
    }

    /**
     * Specifies other arguments which cannot be used in conjunction with tihs argument.  Comma-separated list.
     * @param source Original field specifying command-line arguments.
     * @return A comma-separated list of exclusive arguments, or null if none are present.
     */
    protected String getExclusiveOf( ArgumentSource source ) {
        Argument description = getArgumentDescription(source);
        return description.exclusiveOf().trim().length() > 0 ? description.exclusiveOf().trim() : null;
    }

    /**
     * A regular expression which can be used for validation.
     * @param source Original field specifying command-line arguments.
     * @return a JVM regex-compatible regular expression, or null to permit any possible value.
     */
    protected String getValidationRegex( ArgumentSource source ) {
        Argument description = getArgumentDescription(source);
        return description.validation().trim().length() > 0 ? description.validation().trim() : null;
    }

    /**
     * If the argument source only accepts a small set of options, populate the returned list with
     * those options.  Otherwise, leave the list empty.
     * @param source 
     * @return
     */
    protected List<String> getValidOptions( ArgumentSource source ) {
        if(!source.field.getType().isEnum())
            return null;
        List<String> validOptions = new ArrayList<String>();
        for(Object constant: source.field.getType().getEnumConstants())
            validOptions.add(constant.toString());
        return validOptions;
    }

    /**
     * Gets the value of an argument with the given full name, from the collection of ArgumentMatches.
     * If the argument matches multiple values, an exception will be thrown.
     * @param definition Definition of the argument for which to find matches.
     * @param matches The matches for the given argument.
     * @return The value of the argument if available, or null if not present.
     */
    protected String getArgumentValue( ArgumentDefinition definition, ArgumentMatches matches ) {
        Collection<String> argumentValues = getArgumentValues( definition, matches );
        if( argumentValues.size() > 1 )
            throw new StingException("Multiple values associated with given definition, but this argument expects only one: " + definition.fullName);
        return argumentValues.size() > 0 ? argumentValues.iterator().next() : null;
    }

    /**
     * Gets the values of an argument with the given full name, from the collection of ArgumentMatches.
     * @param definition Definition of the argument for which to find matches.
     * @param matches The matches for the given argument.
     * @return The value of the argument if available, or an empty collection if not present.
     */
    protected Collection<String> getArgumentValues( ArgumentDefinition definition, ArgumentMatches matches ) {
        Collection<String> values = new ArrayList<String>();
        for( ArgumentMatch match: matches ) {
            if( match.definition.equals(definition) )
                values.addAll(match.values());
        }
        return values;
    }

    /**
     * Retrieves the argument description from the given argument source.  Will throw an exception if
     * the given ArgumentSource
     * @param source source of the argument.
     * @return Argument description annotation associated with the given field.
     */
    protected Argument getArgumentDescription( ArgumentSource source ) {
        if( !source.field.isAnnotationPresent(Argument.class) )
            throw new StingException("ArgumentAnnotation is not present for the argument field: " + source.field.getName());
        return source.field.getAnnotation(Argument.class);
    }    
}

/**
 * Parse simple argument types: java primitives, wrapper classes, and anything that has
 * a simple String constructor.
 */
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
    protected Object parse( ArgumentSource source, Class type, ArgumentMatches matches ) {
        String value = getArgumentValue( createDefaultArgumentDefinition(source), matches );

        // lets go through the types we support
        try {
            if (type.isPrimitive()) {
                Method valueOf = primitiveToWrapperMap.get(type).getMethod("valueOf",String.class);
                return valueOf.invoke(null,value.trim());
            } else if (type.isEnum()) {
                Object[] vals = type.getEnumConstants();
                Object defaultEnumeration = null;  // as we look at options, record the default option if it exists
                for (Object val : vals) {
                    if (String.valueOf(val).equalsIgnoreCase(value)) return val;
                    try { if (type.getField(val.toString()).isAnnotationPresent(EnumerationArgumentDefault.class)) defaultEnumeration = val; }
                    catch (NoSuchFieldException e) { throw new StingException("parsing " + type.toString() + "doesn't contain the field " + val.toString()); }
                }
                // if their argument has no value (null), and there's a default, return that default for the enum value
                if (defaultEnumeration != null && value == null)
                    return defaultEnumeration;
                // if their argument has no value and there's no default, throw a missing argument value exception.
                // TODO: Clean this up so that null values never make it to this point.  To fix this, we'll have to clean up the implementation of -U.
                if (value == null)
                    throw new MissingArgumentValueException(Collections.singleton(createDefaultArgumentDefinition(source)));
                throw new UnknownEnumeratedValueException(createDefaultArgumentDefinition(source),value);
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
            throw new StingException("constructFromString:InvocationTargetException: Failed conversion - this is most likely caused by using an incorrect data type (e.g. a double when an int is required)");
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

/**
 * Process compound argument types: arrays, and typed and untyped collections.
 */
class CompoundArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    @Override
    public boolean supports( Class type ) {
        return ( Collection.class.isAssignableFrom(type) || type.isArray() );
    }
    
    @Override
    public Object parse( ArgumentSource source, Class type, ArgumentMatches matches )
    {
        Class componentType = null;

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
            if( source.field.getGenericType() instanceof ParameterizedType) {
                ParameterizedType parameterizedType = (ParameterizedType)source.field.getGenericType();
                if( parameterizedType.getActualTypeArguments().length > 1 )
                    throw new IllegalArgumentException("Unable to determine collection type of field: " + source.field.toString());
                componentType = (Class)parameterizedType.getActualTypeArguments()[0];
            }
            else
                componentType = String.class;

            ArgumentTypeDescriptor componentArgumentParser = ArgumentTypeDescriptor.create( componentType );

            Collection collection = null;
            try {
                collection = (Collection)type.newInstance();
            }
            catch (InstantiationException e) {
                logger.fatal("ArgumentParser: InstantiationException: cannot convert field " + source.field.getName());
                throw new StingException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
            }
            catch (IllegalAccessException e) {
                logger.fatal("ArgumentParser: IllegalAccessException: cannot convert field " + source.field.getName());
                throw new StingException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
            }

            for( ArgumentMatch match: matches ) {
                for( ArgumentMatch value: match )
                    collection.add( componentArgumentParser.parse(source,componentType,new ArgumentMatches(value)) );
            }

            return collection;

        }
        else if( type.isArray() ) {
            componentType = type.getComponentType();
            ArgumentTypeDescriptor componentArgumentParser = ArgumentTypeDescriptor.create( componentType );

            // Assemble a collection of individual values used in this computation.
            Collection<ArgumentMatch> values = new ArrayList<ArgumentMatch>();
            for( ArgumentMatch match: matches ) {
                for( ArgumentMatch value: match )
                    values.add(value);
            }

            Object arr = Array.newInstance(componentType,values.size());

            int i = 0;
            for( ArgumentMatch value: values )
                Array.set( arr,i++,componentArgumentParser.parse(source,componentType,new ArgumentMatches(value)));

            return arr;
        }
        else
            throw new StingException("Unsupported compound argument type: " + type);
    }
}
