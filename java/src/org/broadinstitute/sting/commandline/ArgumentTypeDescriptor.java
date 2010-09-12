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

import org.broadinstitute.sting.utils.exceptions.GATKException;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.gatk.walkers.Multiplex;
import org.broadinstitute.sting.gatk.walkers.Multiplexer;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.lang.annotation.Annotation;
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
    private static Class[] ARGUMENT_ANNOTATIONS = {Input.class, Output.class, Argument.class};

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(ArgumentTypeDescriptor.class);
    
    /**
     * Class reference to the different types of descriptors that the create method can create.
     * The type of set used must be ordered (but not necessarily sorted).
     */
    private static Set<ArgumentTypeDescriptor> descriptors = new LinkedHashSet<ArgumentTypeDescriptor>( Arrays.asList(new SimpleArgumentTypeDescriptor(),
                                                                                                                      new CompoundArgumentTypeDescriptor(),
                                                                                                                      new MultiplexArgumentTypeDescriptor()) );

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
        throw new GATKException("Can't process command-line arguments of type: " + type.getName());
    }

    /**
     * Does this descriptor support classes of the given type?
     * @param type The type to check.
     * @return true if this descriptor supports the given type, false otherwise.
     */
    public abstract boolean supports( Class type );

    /**
     * Returns false if a type-specific default can be employed.
     * @param source Source of the command-line argument.
     * @return True to throw in a type specific default.  False otherwise.
     */
    public boolean createsTypeDefault(ArgumentSource source,Class type) { return false; }

    /**
     * Generates a default for the given type.
     * @param source Source of the command-line argument.
     * @return A default value for the given type.
     */
    public Object createTypeDefault(ArgumentSource source,Class type) { throw new UnsupportedOperationException("Unable to create default for type " + getClass()); }

    /**
     * Given the given argument source and attributes, synthesize argument definitions for command-line arguments.
     * @param source Source class and field for the given argument.
     * @return A list of command-line argument definitions supporting this field.
     */
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source ) {
        return Collections.singletonList(createDefaultArgumentDefinition(source));
    }

    /**
     * Parses an argument source to an object.
     * WARNING!  Mandatory side effect of parsing!  Each parse routine should register the tags it finds with the proper CommandLineProgram.
     * TODO: Fix this, perhaps with an event model indicating that a new argument has been created.
     *
     * @param parsingEngine The engine responsible for parsing.
     * @param source The source used to find the matches.
     * @param matches The matches for the source.
     * @return The parsed object.
     */
    public Object parse(ParsingEngine parsingEngine, ArgumentSource source, ArgumentMatches matches) {
        return parse(parsingEngine, source, source.field.getType(), matches);
    }

    /**
     * Returns true if the field is a collection or an array.
     * @param source The argument source to check.
     * @return true if the field is a collection or an array.
     */
    public boolean isMultiValued( ArgumentSource source ) {
        Class argumentType = source.field.getType();
        return Collection.class.isAssignableFrom(argumentType) || argumentType.isArray();
    }

    /**
     * By default, argument sources create argument definitions with a set of default values.
     * Use this method to create the one simple argument definition.
     * @param source argument source for which to create a default definition.
     * @return The default definition for this argument source.
     */
    protected ArgumentDefinition createDefaultArgumentDefinition( ArgumentSource source ) {
        Annotation argumentAnnotation = getArgumentAnnotation(source);
        return new ArgumentDefinition( ArgumentIOType.getIOType(argumentAnnotation),
                                       source.field.getType(),
                                       ArgumentDefinition.getFullName(argumentAnnotation, source.field.getName()),
                                       ArgumentDefinition.getShortName(argumentAnnotation),
                                       ArgumentDefinition.getDoc(argumentAnnotation),
                                       source.isRequired() && !source.overridesDefault() && !source.isFlag() && !source.isDeprecated(),
                                       source.isFlag(),
                                       source.isMultiValued(),
                                       source.isHidden(),
                                       getCollectionComponentType(source.field),
                                       ArgumentDefinition.getExclusiveOf(argumentAnnotation),
                                       ArgumentDefinition.getValidationRegex(argumentAnnotation),
                                       getValidOptions(source) );
    }

    /**
     * Return the component type of a field, or String.class if the type cannot be found.
     * @param field The reflected field to inspect.
     * @return The parameterized component type, or String.class if the parameterized type could not be found.
     * @throws IllegalArgumentException If more than one parameterized type is found on the field.
     */
    protected Class getCollectionComponentType( Field field ) {
        return null;
    }

    /**
     * Parses the argument matches for a class type into an object.
     * @param source The original argument source used to find the matches.
     * @param type The current class type being inspected.  May not match the argument source.field.getType() if this as a collection for example.
     * @param matches The argument matches for the argument source, or the individual argument match for a scalar if this is being called to help parse a collection.
     * @return The individual parsed object matching the argument match with Class type.
     */
    public abstract Object parse( ParsingEngine parsingEngine, ArgumentSource source, Class type, ArgumentMatches matches );

    /**
     * If the argument source only accepts a small set of options, populate the returned list with
     * those options.  Otherwise, leave the list empty.
     * @param source Original field specifying command-line arguments.
     * @return A list of valid options.
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
     * Returns true if the argument with the given full name exists in the collection of ArgumentMatches.
     * @param definition Definition of the argument for which to find matches.
     * @param matches The matches for the given argument.
     * @return true if the argument is present, or false if not present.
     */
    protected boolean argumentIsPresent( ArgumentDefinition definition, ArgumentMatches matches ) {
        for( ArgumentMatch match: matches ) {
            if( match.definition.equals(definition) )
                return true;
        }
        return false;
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
            throw new UserException.CommandLineException("Multiple values associated with given definition, but this argument expects only one: " + definition.fullName);
        return argumentValues.size() > 0 ? argumentValues.iterator().next() : null;
    }

    /**
     * Gets the tags associated with a given command-line argument.
     * If the argument matches multiple values, an exception will be thrown.
     * @param matches The matches for the given argument.
     * @return The value of the argument if available, or null if not present.
     */
    protected List<String> getArgumentTags(ArgumentMatches matches) {
        List<String> tags = new ArrayList<String>();
        for( ArgumentMatch match: matches ) {
                tags.addAll(match.tags);
        }
        return tags;
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
    @SuppressWarnings("unchecked")
    protected static Annotation getArgumentAnnotation( ArgumentSource source ) {
        for (Class annotation: ARGUMENT_ANNOTATIONS)
            if (source.field.isAnnotationPresent(annotation))
                return source.field.getAnnotation(annotation);
        throw new GATKException("ArgumentAnnotation is not present for the argument field: " + source.field.getName());
    }

    /**
     * Returns true if an argument annotation is present
     * @param field The field to check for an annotation.
     * @return True if an argument annotation is present on the field.
     */
    @SuppressWarnings("unchecked")
    public static boolean isArgumentAnnotationPresent(Field field) {
        for (Class annotation: ARGUMENT_ANNOTATIONS)
            if (field.isAnnotationPresent(annotation))
                return true;
        return false;
    }

    /**
     * Returns true if the given annotation is hidden from the help system.
     * @param field Field to test.
     * @return True if argument should be hidden.  False otherwise.
     */
    public static boolean isArgumentHidden(Field field) {
        return field.isAnnotationPresent(Hidden.class);
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
    public Object parse(ParsingEngine parsingEngine, ArgumentSource source, Class type, ArgumentMatches matches) {
        if (source.isFlag())
            return true;
        ArgumentDefinition defaultDefinition = createDefaultArgumentDefinition(source);
        String value = getArgumentValue( defaultDefinition, matches );
        Object result;
        List<String> tags = getArgumentTags( matches );

        // lets go through the types we support
        try {
            if (type.isPrimitive()) {
                Method valueOf = primitiveToWrapperMap.get(type).getMethod("valueOf",String.class);
                result = valueOf.invoke(null,value.trim());
            } else if (type.isEnum()) {
                Object[] vals = type.getEnumConstants();
                Object defaultEnumeration = null;  // as we look at options, record the default option if it exists
                for (Object val : vals) {
                    if (String.valueOf(val).equalsIgnoreCase(value)) return val;
                    try { if (type.getField(val.toString()).isAnnotationPresent(EnumerationArgumentDefault.class)) defaultEnumeration = val; }
                    catch (NoSuchFieldException e) { throw new GATKException("parsing " + type.toString() + "doesn't contain the field " + val.toString()); }
                }
                // if their argument has no value (null), and there's a default, return that default for the enum value
                if (defaultEnumeration != null && value == null)
                    result = defaultEnumeration;
                // if their argument has no value and there's no default, throw a missing argument value exception.
                // TODO: Clean this up so that null values never make it to this point.  To fix this, we'll have to clean up the implementation of -U.
                else if (value == null)
                    throw new MissingArgumentValueException(Collections.singleton(createDefaultArgumentDefinition(source)));
                else
                    throw new UnknownEnumeratedValueException(createDefaultArgumentDefinition(source),value);
            } else {
                Constructor ctor = type.getConstructor(String.class);
                result = ctor.newInstance(value);
            }
        } catch (InvocationTargetException e) {
            throw new GATKException("constructFromString:InvocationTargetException: Failed conversion - this is most likely caused by using an incorrect data type (e.g. a double when an int is required)");
        } catch (Exception e) {
            throw new DynamicClassResolutionException(String.class, e);
        }


        // WARNING: Side effect!
        parsingEngine.addTags(result,tags);

        return result;
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
    @SuppressWarnings("unchecked")
    public Object parse(ParsingEngine parsingEngine,ArgumentSource source, Class type, ArgumentMatches matches) {
        Class componentType;
        Object result;
        List<String> tags = new ArrayList<String>();

        if( Collection.class.isAssignableFrom(type) ) {

            // If this is a generic interface, pick a concrete implementation to create and pass back.
            // Because of type erasure, don't worry about creating one of exactly the correct type.
            if( Modifier.isInterface(type.getModifiers()) || Modifier.isAbstract(type.getModifiers()) )
            {
                if( java.util.List.class.isAssignableFrom(type) ) type = ArrayList.class;
                else if( java.util.Queue.class.isAssignableFrom(type) ) type = java.util.ArrayDeque.class;
                else if( java.util.Set.class.isAssignableFrom(type) ) type = java.util.TreeSet.class;
            }

            componentType = getCollectionComponentType( source.field );
            ArgumentTypeDescriptor componentArgumentParser = ArgumentTypeDescriptor.create( componentType );

            Collection collection;
            try {
                collection = (Collection)type.newInstance();
            }
            catch (InstantiationException e) {
                logger.fatal("ArgumentParser: InstantiationException: cannot convert field " + source.field.getName());
                throw new GATKException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
            }
            catch (IllegalAccessException e) {
                logger.fatal("ArgumentParser: IllegalAccessException: cannot convert field " + source.field.getName());
                throw new GATKException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
            }

            for( ArgumentMatch match: matches ) {
                for( ArgumentMatch value: match ) {
                    collection.add( componentArgumentParser.parse(parsingEngine,source,componentType,new ArgumentMatches(value)) );
                    tags.addAll(value.tags);
                }
            }

            result = collection;

        }
        else if( type.isArray() ) {
            componentType = type.getComponentType();
            ArgumentTypeDescriptor componentArgumentParser = ArgumentTypeDescriptor.create( componentType );

            // Assemble a collection of individual values used in this computation.
            Collection<ArgumentMatch> values = new ArrayList<ArgumentMatch>();
            for( ArgumentMatch match: matches )
                for( ArgumentMatch value: match )
                    values.add(value);

            result = Array.newInstance(componentType,values.size());

            int i = 0;
            for( ArgumentMatch value: values ) {
                Array.set( result,i++,componentArgumentParser.parse(parsingEngine,source,componentType,new ArgumentMatches(value)));
                tags.addAll(value.tags);
            }
        }
        else
            throw new GATKException("Unsupported compound argument type: " + type);

        // WARNING: Side effect!
        parsingEngine.addTags(result,tags);

        return result;
    }

    /**
     * Return the component type of a field, or String.class if the type cannot be found.
     * @param field The reflected field to inspect.
     * @return The parameterized component type, or String.class if the parameterized type could not be found.
     * @throws IllegalArgumentException If more than one parameterized type is found on the field.
     */
    @Override
    protected Class getCollectionComponentType( Field field ) {
            // If this is a parameterized collection, find the contained type.  If blow up if more than one type exists.
            if( field.getGenericType() instanceof ParameterizedType) {
                ParameterizedType parameterizedType = (ParameterizedType)field.getGenericType();
                if( parameterizedType.getActualTypeArguments().length > 1 )
                    throw new IllegalArgumentException("Unable to determine collection type of field: " + field.toString());
                return (Class)parameterizedType.getActualTypeArguments()[0];
            }
            else
                return String.class;
    }
}

class MultiplexArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * The multiplexer controlling how data is split.
     */
    private final Multiplexer multiplexer;

    /**
     * The set of identifiers for the multiplexed entries.
     */
    private final Collection<?> multiplexedIds;

    public MultiplexArgumentTypeDescriptor() {
        this.multiplexer = null;
        this.multiplexedIds = null;
    }

    /**
     * Private constructor to use in creating a closure of the MultiplexArgumentTypeDescriptor specific to the
     * given set of multiplexed ids.
     * @param multiplexedIds The collection of multiplexed entries
     */
    private MultiplexArgumentTypeDescriptor(final Multiplexer multiplexer, final Collection<?> multiplexedIds) {
        this.multiplexer = multiplexer;
        this.multiplexedIds = multiplexedIds;
    }

    @Override
    public boolean supports( Class type ) {
        return ( Map.class.isAssignableFrom(type) );
    }

    @Override
    public boolean createsTypeDefault(ArgumentSource source,Class type) {
        // Multiplexing always creates a type default.
        return true;
    }

    @Override
    public Object createTypeDefault(ArgumentSource source,Class type) {
        if(multiplexer == null || multiplexedIds == null)
            throw new GATKException("No multiplexed ids available");

        Map<Object,Object> multiplexedMapping = new HashMap<Object,Object>();
        Class componentType = getCollectionComponentType(source.field);
        ArgumentTypeDescriptor componentTypeDescriptor = ArgumentTypeDescriptor.create(componentType);

        for(Object id: multiplexedIds) {
            Object value = null;
            if(componentTypeDescriptor.createsTypeDefault(source,componentType))
                value = componentTypeDescriptor.createTypeDefault(source,componentType);
            multiplexedMapping.put(id,value);
        }
        return multiplexedMapping;
    }


    @Override
    public Object parse(ParsingEngine parsingEngine, ArgumentSource source, Class type, ArgumentMatches matches) {
        if(multiplexedIds == null)
            throw new GATKException("Cannot directly parse a MultiplexArgumentTypeDescriptor; must create a derivative type descriptor first.");

        Map<Object,Object> multiplexedMapping = new HashMap<Object,Object>();

        Class componentType = getCollectionComponentType(source.field);


        for(Object id: multiplexedIds) {
            Object value = ArgumentTypeDescriptor.create(componentType).parse(parsingEngine,source,componentType,matches.transform(multiplexer,id));
            multiplexedMapping.put(id,value);
        }

        parsingEngine.addTags(multiplexedMapping,getArgumentTags(matches));

        return multiplexedMapping;
    }

    public MultiplexArgumentTypeDescriptor createCustomTypeDescriptor(ArgumentSource dependentArgument,Object containingObject) {
        String[] sourceFields = dependentArgument.field.getAnnotation(Multiplex.class).arguments();

        List<ArgumentSource> allSources = ParsingEngine.extractArgumentSources(containingObject.getClass());
        Class[] sourceTypes = new Class[sourceFields.length];
        Object[] sourceValues = new Object[sourceFields.length];
        int currentField = 0;

        for(String sourceField: sourceFields) {
            boolean fieldFound = false;
            for(ArgumentSource source: allSources) {
                if(!source.field.getName().equals(sourceField))
                    continue;
                if(source.field.isAnnotationPresent(Multiplex.class))
                    throw new GATKException("Command-line arguments can only depend on independent fields");
                sourceTypes[currentField] = source.field.getType();
                sourceValues[currentField] = JVMUtils.getFieldValue(source.field,containingObject);
                currentField++;
                fieldFound = true;
            }
            if(!fieldFound)
                throw new GATKException(String.format("Unable to find source field %s, referred to by dependent field %s",sourceField,dependentArgument.field.getName()));
        }

        Class<? extends Multiplexer> multiplexerType = dependentArgument.field.getAnnotation(Multiplex.class).value();
        Constructor<? extends Multiplexer> multiplexerConstructor = null;
        try {
            multiplexerConstructor = multiplexerType.getConstructor(sourceTypes);
            multiplexerConstructor.setAccessible(true);
        }
        catch(NoSuchMethodException ex) {
            throw new GATKException(String.format("Unable to find constructor for class %s with parameters %s",multiplexerType.getName(),Arrays.deepToString(sourceFields)),ex);
        }

        Multiplexer multiplexer = null;
        try {
            multiplexer = multiplexerConstructor.newInstance(sourceValues);
        }
        catch(IllegalAccessException ex) {
            throw new GATKException(String.format("Constructor for class %s with parameters %s is inaccessible",multiplexerType.getName(),Arrays.deepToString(sourceFields)),ex);
        }
        catch(InstantiationException ex) {
            throw new GATKException(String.format("Can't create class %s with parameters %s",multiplexerType.getName(),Arrays.deepToString(sourceFields)),ex);
        }
        catch(InvocationTargetException ex) {
            throw new GATKException(String.format("Can't invoke constructor of class %s with parameters %s",multiplexerType.getName(),Arrays.deepToString(sourceFields)),ex);
        }

        return new MultiplexArgumentTypeDescriptor(multiplexer,multiplexer.multiplex());
    }

    /**
     * Return the component type of a field, or String.class if the type cannot be found.
     * @param field The reflected field to inspect.
     * @return The parameterized component type, or String.class if the parameterized type could not be found.
     * @throws IllegalArgumentException If more than one parameterized type is found on the field.
     */
    @Override
    protected Class getCollectionComponentType( Field field ) {
        // Multiplex arguments must resolve to maps from which the clp should extract the second type.
        if( field.getGenericType() instanceof ParameterizedType) {
            ParameterizedType parameterizedType = (ParameterizedType)field.getGenericType();
            if( parameterizedType.getActualTypeArguments().length != 2 )
                throw new IllegalArgumentException("Unable to determine collection type of field: " + field.toString());
            return (Class)parameterizedType.getActualTypeArguments()[1];
        }
        else
            return String.class;
    }
}
