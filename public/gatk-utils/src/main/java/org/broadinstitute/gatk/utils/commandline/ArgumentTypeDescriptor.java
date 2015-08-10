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

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.classloader.JVMUtils;
import org.broadinstitute.gatk.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.refdata.tracks.FeatureManager;
import org.broadinstitute.gatk.utils.text.XReadLines;

import java.io.File;
import java.io.IOException;
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
     * our log, which we want to capture anything from org.broadinstitute.gatk
     */
    protected static final Logger logger = Logger.getLogger(ArgumentTypeDescriptor.class);

    /**
     * Fetch the given descriptor from the descriptor repository.
     * @param descriptors the descriptors from which to select a good match.
     * @param type Class for which to specify a descriptor.
     * @return descriptor for the given type.
     */
    public static ArgumentTypeDescriptor selectBest( Collection<ArgumentTypeDescriptor> descriptors, Class type ) {
        for( ArgumentTypeDescriptor descriptor: descriptors ) {
            if( descriptor.supports(type) )
                return descriptor;
        }
        throw new ReviewedGATKException("Can't process command-line arguments of type: " + type.getName());
    }

    /**
     * Returns true if the file will be compressed.
     * @param writerFileName Name of the file
     * @return true if the file will be compressed.
     */
    public static boolean isCompressed(String writerFileName) {
        return writerFileName != null && AbstractFeatureReader.hasBlockCompressedExtension(writerFileName);
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
    public boolean createsTypeDefault(ArgumentSource source) { return false; }

    /**
     * Returns a documentation-friendly value for the default of a type descriptor.
     * Must be overridden if createsTypeDefault return true.  cannot be called otherwise
     * @param source Source of the command-line argument.
     * @return Friendly string of the default value, for documentation.  If doesn't create a default, throws
     * and UnsupportedOperationException
     */
    public String typeDefaultDocString(ArgumentSource source) {
        throw new UnsupportedOperationException();
    }

    /**
     * Generates a default for the given type.
     *
     * @param parsingEngine the parsing engine used to validate this argument type descriptor.
     * @param source Source of the command-line argument.
     * @param type Type of value to create, in case the command-line argument system wants influence.
     * @return A default value for the given type.
     */
    public Object createTypeDefault(ParsingEngine parsingEngine,ArgumentSource source, Type type) { throw new UnsupportedOperationException("Unable to create default for type " + getClass()); }

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
        return parse(parsingEngine, source, source.field.getGenericType(), matches);
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
                source.isRequired() && !createsTypeDefault(source) && !source.isFlag() && !source.isDeprecated(),
                source.isFlag(),
                source.isMultiValued(),
                source.isHidden(),
                makeRawTypeIfNecessary(getCollectionComponentType(source.field)),
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
    protected Type getCollectionComponentType( Field field ) {
        return null;
    }

    /**
     * Parses the argument matches for a class type into an object.
     * @param source The original argument source used to find the matches.
     * @param type The current class type being inspected.  May not match the argument source.field.getType() if this as a collection for example.
     * @param matches The argument matches for the argument source, or the individual argument match for a scalar if this is being called to help parse a collection.
     * @return The individual parsed object matching the argument match with Class type.
     */
    public abstract Object parse( ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches );

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
    protected ArgumentMatchValue getArgumentValue( ArgumentDefinition definition, ArgumentMatches matches ) {
        Collection<ArgumentMatchValue> argumentValues = getArgumentValues( definition, matches );
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
    protected Tags getArgumentTags(ArgumentMatches matches) {
        Tags tags = new Tags();
        for(ArgumentMatch match: matches) {
            if(!tags.isEmpty() && !match.tags.isEmpty())
                throw new ReviewedGATKException("BUG: multiple conflicting sets of tags are available, and the type descriptor specifies no way of resolving the conflict.");
            tags = match.tags;
        }
        return tags;
    }

    /**
     * Gets the values of an argument with the given full name, from the collection of ArgumentMatches.
     * @param definition Definition of the argument for which to find matches.
     * @param matches The matches for the given argument.
     * @return The value of the argument if available, or an empty collection if not present.
     */
    protected Collection<ArgumentMatchValue> getArgumentValues( ArgumentDefinition definition, ArgumentMatches matches ) {
        Collection<ArgumentMatchValue> values = new ArrayList<ArgumentMatchValue>();
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
        throw new ReviewedGATKException("ArgumentAnnotation is not present for the argument field: " + source.field.getName());
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

    public static Class makeRawTypeIfNecessary(Type t) {
        if ( t == null )
            return null;
        else if ( t instanceof ParameterizedType )
            return (Class)((ParameterizedType) t).getRawType();
        else if ( t instanceof Class ) {
            return (Class)t;
        } else {
            throw new IllegalArgumentException("Unable to determine Class-derived component type of field: " + t);
        }
    }

    /**
     * The actual argument parsing method.
     * @param source             source
     * @param type               type to check
     * @param matches            matches
     * @param tags               argument tags
     * @return the RodBinding/IntervalBinding object depending on the value of createIntervalBinding.
     */
    protected Object parseBinding(ArgumentSource source, Type type, ArgumentMatches matches, Tags tags) {
        ArgumentDefinition defaultDefinition = createDefaultArgumentDefinition(source);
        ArgumentMatchValue value = getArgumentValue(defaultDefinition, matches);
        @SuppressWarnings("unchecked")
        Class<? extends Feature> parameterType = JVMUtils.getParameterizedTypeClass(type);
        String name = defaultDefinition.fullName;

        return parseBinding(value, parameterType, type, name, tags, source.field.getName());
    }

    /**
     *
     * @param value The source of the binding
     * @param parameterType The Tribble Feature parameter type
     * @param bindingClass The class type for the binding (ex: RodBinding, IntervalBinding, etc.) Must have the correct constructor for creating the binding.
     * @param bindingName The name of the binding passed to the constructor.
     * @param tags Tags for the binding used for parsing and passed to the constructor.
     * @param fieldName The name of the field that was parsed. Used for error reporting.
     * @return The newly created binding object of type bindingClass.
     */
    public static Object parseBinding(ArgumentMatchValue value, Class<? extends Feature> parameterType, Type bindingClass,
                                      String bindingName, Tags tags, String fieldName) {
        try {
            String tribbleType = null;
            // must have one or two tag values here
            if ( tags.getPositionalTags().size() > 2 ) {
                throw new UserException.CommandLineException(
                        String.format("Unexpected number of positional tags for argument %s : %s. " +
                                "Rod bindings only support -X:type and -X:name,type argument styles",
                                value.asString(), fieldName));
            } else if ( tags.getPositionalTags().size() == 2 ) {
                // -X:name,type style
                bindingName = tags.getPositionalTags().get(0);
                tribbleType = tags.getPositionalTags().get(1);

                FeatureManager manager = new FeatureManager();
                if ( manager.getByName(tribbleType) == null )
                    throw new UserException.UnknownTribbleType(
                            tribbleType,
                            String.format("Unable to find tribble type '%s' provided on the command line. " +
                                    "Please select a correct type from among the supported types:%n%s",
                                    tribbleType, manager.userFriendlyListOfAvailableFeatures(parameterType)));

            } else {
                // case with 0 or 1 positional tags
                FeatureManager manager = new FeatureManager();

                // -X:type style is a type when we cannot determine the type dynamically
                String tag1 = tags.getPositionalTags().size() == 1 ? tags.getPositionalTags().get(0) : null;
                if ( tag1 != null ) {
                    if ( manager.getByName(tag1) != null ) // this a type
                        tribbleType = tag1;
                    else
                        bindingName = tag1;
                }

                if ( tribbleType == null ) {
                    // try to determine the file type dynamically
                    File file = value.asFile();
                    if ( file.canRead() && file.isFile() ) {
                        FeatureManager.FeatureDescriptor featureDescriptor = manager.getByFiletype(file);
                        if ( featureDescriptor != null ) {
                            tribbleType = featureDescriptor.getName();
                            logger.debug("Dynamically determined type of " + file + " to be " + tribbleType);
                        }
                    }

                    if ( tribbleType == null ) {
                        // IntervalBinding can be created from a normal String
                        Class rawType = (makeRawTypeIfNecessary(bindingClass));
                        try {
                            return rawType.getConstructor(String.class).newInstance(value.asString());
                        } catch (NoSuchMethodException e) {
                            /* ignore */
                        }

                        if ( ! file.exists() ) {
                            throw new UserException.CouldNotReadInputFile(file, "file \'"+ file + "\' does not exist");
                        } else if ( ! file.canRead() || ! file.isFile() ) {
                            throw new UserException.CouldNotReadInputFile(file, "file \'"+ file + "\' could not be read");
                        } else {
                            throw new UserException.CommandLineException(
                                    String.format("No tribble type was provided on the command line and the type of the file \'"+ file + "\' could not be determined dynamically. " +
                                                    "Please add an explicit type tag :NAME listing the correct type from among the supported types:%n%s",
                                            manager.userFriendlyListOfAvailableFeatures(parameterType)));
                        }
                    }
                }
            }

            Constructor ctor = (makeRawTypeIfNecessary(bindingClass)).getConstructor(Class.class, String.class, String.class, String.class, Tags.class);
            return ctor.newInstance(parameterType, bindingName, value.asString(), tribbleType, tags);
        } catch (Exception e) {
            if ( e instanceof UserException )
                throw ((UserException)e);
            else
                throw new UserException.CommandLineException(
                        String.format("Failed to parse value %s for argument %s. Message: %s",
                                value, fieldName, e.getMessage()));
        }
    }

    /**
     * Parse the source of a RodBindingCollection, which can be either a file of RodBindings or an actual RodBinding.
     *
     * @param parsingEngine the parsing engine used to validate this argument type descriptor
     * @param source             source
     * @param type               type
     * @param matches            matches
     * @param tags               argument tags
     * @return the newly created binding object
     */
    public Object parseRodBindingCollectionSource(final ParsingEngine parsingEngine,
                                                  final ArgumentSource source,
                                                  final Type type,
                                                  final ArgumentMatches matches,
                                                  final Tags tags) {

        final ArgumentDefinition defaultDefinition = createDefaultArgumentDefinition(source);
        final ArgumentMatchValue value = getArgumentValue(defaultDefinition, matches);
        @SuppressWarnings("unchecked")
        Class<? extends Feature> parameterType = JVMUtils.getParameterizedTypeClass(type);
        String name = defaultDefinition.fullName;

        // if this a list of files, get those bindings
        final File file = value.asFile();
        try {
            if (file.getAbsolutePath().endsWith(".list")) {
                return getRodBindingsCollection(file, parsingEngine, parameterType, name, tags, source.field.getName());
            }
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }

        // otherwise, treat this as an individual binding
        final RodBinding binding = (RodBinding)parseBinding(value, parameterType, RodBinding.class, name, tags, source.field.getName());
        parsingEngine.addTags(binding, tags);
        parsingEngine.addRodBinding(binding);
        return RodBindingCollection.createRodBindingCollectionOfType(parameterType, Arrays.asList(binding));
    }

    /**
     * Retrieve and parse a collection of RodBindings from the given file.
     *
     * If the file contains duplicate entries or is empty, an exception will be thrown.
     *
     * @param file             the source file
     * @param parsingEngine    the engine responsible for parsing
     * @param parameterType    the Tribble Feature parameter type
     * @param bindingName      the name of the binding passed to the constructor.
     * @param defaultTags      general tags for the binding used for parsing and passed to the constructor.
     * @param fieldName        the name of the field that was parsed. Used for error reporting.
     * @return the newly created collection of binding objects.
     */
    public static Object getRodBindingsCollection(final File file,
                                                  final ParsingEngine parsingEngine,
                                                  final Class<? extends Feature> parameterType,
                                                  final String bindingName,
                                                  final Tags defaultTags,
                                                  final String fieldName) throws IOException {
        final List<RodBinding> bindings = new ArrayList<>();

        // Keep track of the files in this list so that we can check for duplicates and empty files
        final Set<String> fileValues = new HashSet<>();

        // parse each line separately using the given Tags if none are provided on each line
        for ( final String line: new XReadLines(file) ) {
            final String[] tokens = line.split("\\s+");
            final RodBinding binding;

            if ( tokens.length == 0 ) {
                continue; // empty line, so do nothing
            }
            // use the default tags if none are provided for this binding
            else if ( tokens.length == 1 ) {
                final ArgumentMatchValue value = parseAndValidateArgumentMatchValue(tokens[0], fileValues, fieldName, file.getName());
                binding = (RodBinding)parseBinding(value, parameterType, RodBinding.class, bindingName, defaultTags, fieldName);
                parsingEngine.addTags(binding, defaultTags);

            }
            // use the new tags if provided
            else if ( tokens.length == 2 ) {
                final Tags tags = ParsingMethod.parseTags(fieldName, tokens[0]);
                final ArgumentMatchValue value = parseAndValidateArgumentMatchValue(tokens[1], fileValues, fieldName, file.getName());
                binding = (RodBinding)parseBinding(value, parameterType, RodBinding.class, bindingName, tags, fieldName);
                parsingEngine.addTags(binding, tags);
            } else {
                throw new UserException.BadArgumentValue(fieldName, "data lines should consist of an optional set of tags along with a path to a file; too many tokens are present for line: " + line);
            }

            bindings.add(binding);
            parsingEngine.addRodBinding(binding);
        }

        if (fileValues.isEmpty()) {
            throw new UserException.BadArgumentValue(fieldName, "The input list " + file.getName() + " is empty.");
        }

        return RodBindingCollection.createRodBindingCollectionOfType(parameterType, bindings);
    }

    /**
     * Validates the resource file name and constructs an ArgumentMatchValue from it.
     *
     * If the list name has already been processed in the current list, throws a UserException, otherwise
     * creates an ArgumentMatchValue to represent the list.
     *
     * @param token Name of the ROD resource file.
     * @param fileValues Set of names of ROD files that have already been processed.
     * @param fieldName Name of the argument field being populated.
     * @param listFileName Name of the list file being processed.
     * @return
     */
    private static ArgumentMatchValue parseAndValidateArgumentMatchValue(final String token, final Set<String> fileValues, final String fieldName,
                                                                         final String listFileName) {
        checkForDuplicateFileName(token, fileValues, fieldName, listFileName);
        return new ArgumentMatchStringValue(token);
    }

    /**
     * Checks to make sure that the current file name to be processed has not already been processed.
     *
     * Checks the name of the current file against the names that have already been processed, throwing
     * an informative BadArgumentValue exception if it has already been seen. As a side effect adds the
     * current file name to the set of filenames that have already been processed.
     *
     * @param currentFile Name of the current file to process
     * @param processedFiles Set of file names that have already been processed
     * @param fieldName Name of the argument that is being populated
     * @param listName Filename of the list that is being processed
     */
    protected static void checkForDuplicateFileName(final String currentFile, final Set<String> processedFiles,
                                                    final String fieldName, final String listName) {
        if (processedFiles.contains(currentFile)) {
            throw new UserException.BadArgumentValue(fieldName, "The input list " + listName + " contains file " + currentFile +
                                                     " multiple times, which isn't allowed. If you are intentionally trying to " +
                                                     "include the same file more than once, you will need to specify it in separate file lists.");
        }
        processedFiles.add(currentFile);
    }
}

/**
 * Parser for RodBinding objects
 */
class RodBindingArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * We only want RodBinding class objects
     * @param type The type to check.
     * @return true if the provided class is a RodBinding.class
     */
    @Override
    public boolean supports( Class type ) {
        return isRodBinding(type);
    }

    public static boolean isRodBinding( Class type ) {
        return RodBinding.class.isAssignableFrom(type);
    }

    @Override
    public boolean createsTypeDefault(ArgumentSource source) { return ! source.isRequired(); }

    @Override
    @SuppressWarnings("unchecked")
    public Object createTypeDefault(ParsingEngine parsingEngine, ArgumentSource source, Type type) {
        Class parameterType = JVMUtils.getParameterizedTypeClass(type);
        return RodBinding.makeUnbound((Class<? extends Feature>)parameterType);
    }

    @Override
    public String typeDefaultDocString(ArgumentSource source) {
        return "none";
    }

    @Override
    public Object parse(ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches) {
        Tags tags = getArgumentTags(matches);
        RodBinding rbind = (RodBinding)parseBinding(source, type, matches, tags);
        parsingEngine.addTags(rbind, tags);
        parsingEngine.addRodBinding(rbind);
        return rbind;
    }
}

/**
 * Parser for IntervalBinding objects
 */
class IntervalBindingArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * We only want IntervalBinding class objects
     * @param type The type to check.
     * @return true if the provided class is an IntervalBinding.class
     */
    @Override
    public boolean supports( Class type ) {
        return isIntervalBinding(type);
    }

    public static boolean isIntervalBinding( Class type ) {
        return IntervalBinding.class.isAssignableFrom(type);
    }

    /**
     * See note from RodBindingArgumentTypeDescriptor.parse().
     *
     * @param parsingEngine      parsing engine
     * @param source             source
     * @param type               type to check
     * @param matches            matches
     * @return the IntervalBinding object.
     */
    @Override
    public Object parse(ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches) {
        return parseBinding(source, type, matches, getArgumentTags(matches));
    }
}

/**
 * Parser for RodBindingCollection objects
 */
class RodBindingCollectionArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * We only want RodBindingCollection class objects
     * @param type The type to check.
     * @return true if the provided class is an RodBindingCollection.class
     */
    @Override
    public boolean supports( final Class type ) {
        return isRodBindingCollection(type);
    }

    public static boolean isRodBindingCollection( final Class type ) {
        return RodBindingCollection.class.isAssignableFrom(type);
    }

    /**
     * See note from RodBindingArgumentTypeDescriptor.parse().
     *
     * @param parsingEngine      parsing engine
     * @param source             source
     * @param type               type to check
     * @param matches            matches
     * @return the IntervalBinding object.
     */
    @Override
    public Object parse(final ParsingEngine parsingEngine, final ArgumentSource source, final Type type, final ArgumentMatches matches) {
        final Tags tags = getArgumentTags(matches);
        return parseRodBindingCollectionSource(parsingEngine, source, type, matches, tags);
    }
}

/**
 * Parse simple argument types: java primitives, wrapper classes, and anything that has
 * a simple String constructor.
 */
class SimpleArgumentTypeDescriptor extends ArgumentTypeDescriptor {

    /**
     * @param type  the class type
     * @return true if this class is a binding type, false otherwise
     */
    private boolean isBinding(final Class type) {
        return RodBindingArgumentTypeDescriptor.isRodBinding(type) ||
                IntervalBindingArgumentTypeDescriptor.isIntervalBinding(type) ||
                RodBindingCollectionArgumentTypeDescriptor.isRodBindingCollection(type);
    }


    @Override
    public boolean supports( Class type ) {
        if ( isBinding(type) ) return false;
        if ( type.isPrimitive() ) return true;
        if ( type.isEnum() ) return true;
        if ( primitiveToWrapperMap.containsValue(type) ) return true;

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
    public Object parse(ParsingEngine parsingEngine, ArgumentSource source, Type fulltype, ArgumentMatches matches) {
        Class type = makeRawTypeIfNecessary(fulltype);
        if (source.isFlag())
            return true;

        ArgumentDefinition defaultDefinition = createDefaultArgumentDefinition(source);
        ArgumentMatchValue value = getArgumentValue(defaultDefinition, matches);
        Object result;
        Tags tags = getArgumentTags(matches);

        // lets go through the types we support
        try {
            if (type.isPrimitive()) {
                Method valueOf = primitiveToWrapperMap.get(type).getMethod("valueOf",String.class);
                if(value == null)
                    throw new MissingArgumentValueException(createDefaultArgumentDefinition(source));
                result = valueOf.invoke(null,value.asString().trim());
            } else if (type.isEnum()) {
                Object[] vals = type.getEnumConstants();
                Object defaultEnumeration = null;  // as we look at options, record the default option if it exists
                for (Object val : vals) {
                    if (String.valueOf(val).equalsIgnoreCase(value == null ? null : value.asString())) return val;
                    try { if (type.getField(val.toString()).isAnnotationPresent(EnumerationArgumentDefault.class)) defaultEnumeration = val; }
                    catch (NoSuchFieldException e) { throw new ReviewedGATKException("parsing " + type.toString() + "doesn't contain the field " + val.toString()); }
                }
                // if their argument has no value (null), and there's a default, return that default for the enum value
                if (defaultEnumeration != null && value == null)
                    result = defaultEnumeration;
                    // if their argument has no value and there's no default, throw a missing argument value exception.
                    // TODO: Clean this up so that null values never make it to this point.  To fix this, we'll have to clean up the implementation of -U.
                else if (value == null)
                    throw new MissingArgumentValueException(createDefaultArgumentDefinition(source));
                else
                    throw new UnknownEnumeratedValueException(createDefaultArgumentDefinition(source),value.asString());
            } else if (type.equals(File.class)) {
                result = value == null ? null : value.asFile();
            } else {
                if (value == null)
                    throw new MissingArgumentValueException(createDefaultArgumentDefinition(source));
                Constructor ctor = type.getConstructor(String.class);
                result = ctor.newInstance(value.asString());
            }
        } catch (UserException e) {
            throw e;
        } catch (InvocationTargetException e) {
            throw new UserException.CommandLineException(String.format("Failed to parse value %s for argument %s.  This is most commonly caused by providing an incorrect data type (e.g. a double when an int is required)",
                    value, source.field.getName()));
        } catch (Exception e) {
            throw new DynamicClassResolutionException(String.class, e);
        }

        // TODO FIXME!

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
    public Object parse(ParsingEngine parsingEngine,ArgumentSource source, Type fulltype, ArgumentMatches matches) {
        Class type = makeRawTypeIfNecessary(fulltype);
        Type componentType;
        Object result;

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
            ArgumentTypeDescriptor componentArgumentParser = parsingEngine.selectBestTypeDescriptor(makeRawTypeIfNecessary(componentType));

            Collection collection;
            try {
                collection = (Collection)type.newInstance();
            }
            catch (InstantiationException e) {
                logger.fatal("ArgumentParser: InstantiationException: cannot convert field " + source.field.getName());
                throw new ReviewedGATKException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
            }
            catch (IllegalAccessException e) {
                logger.fatal("ArgumentParser: IllegalAccessException: cannot convert field " + source.field.getName());
                throw new ReviewedGATKException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
            }

            for( ArgumentMatch match: matches ) {
                for( ArgumentMatch value: match ) {
                    Object object = componentArgumentParser.parse(parsingEngine,source,componentType,new ArgumentMatches(value));
                    collection.add( object );
                    // WARNING: Side effect!
                    parsingEngine.addTags(object,value.tags);
                }
            }

            result = collection;

        }
        else if( type.isArray() ) {
            componentType = type.getComponentType();
            ArgumentTypeDescriptor componentArgumentParser = parsingEngine.selectBestTypeDescriptor(makeRawTypeIfNecessary(componentType));

            // Assemble a collection of individual values used in this computation.
            Collection<ArgumentMatch> values = new ArrayList<ArgumentMatch>();
            for( ArgumentMatch match: matches )
                for( ArgumentMatch value: match )
                    values.add(value);

            result = Array.newInstance(makeRawTypeIfNecessary(componentType),values.size());

            int i = 0;
            for( ArgumentMatch value: values ) {
                Object object = componentArgumentParser.parse(parsingEngine,source,componentType,new ArgumentMatches(value));
                Array.set(result,i++,object);
                // WARNING: Side effect!
                parsingEngine.addTags(object,value.tags);
            }
        }
        else
            throw new ReviewedGATKException("Unsupported compound argument type: " + type);

        return result;
    }

    /**
     * Return the component type of a field, or String.class if the type cannot be found.
     * @param field The reflected field to inspect.
     * @return The parameterized component type, or String.class if the parameterized type could not be found.
     * @throws IllegalArgumentException If more than one parameterized type is found on the field.
     */
    @Override
    protected Type getCollectionComponentType( Field field ) {
        // If this is a parameterized collection, find the contained type.  If blow up if more than one type exists.
        if( field.getGenericType() instanceof ParameterizedType) {
            ParameterizedType parameterizedType = (ParameterizedType)field.getGenericType();
            if( parameterizedType.getActualTypeArguments().length > 1 )
                throw new IllegalArgumentException("Unable to determine collection type of field: " + field.toString());
            return parameterizedType.getActualTypeArguments()[0];
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
    public boolean createsTypeDefault(ArgumentSource source) {
        // Multiplexing always creates a type default.
        return true;
    }

    @Override
    public Object createTypeDefault(ParsingEngine parsingEngine,ArgumentSource source, Type type) {
        if(multiplexer == null || multiplexedIds == null)
            throw new ReviewedGATKException("No multiplexed ids available");

        Map<Object,Object> multiplexedMapping = new HashMap<Object,Object>();
        Class componentType = makeRawTypeIfNecessary(getCollectionComponentType(source.field));
        ArgumentTypeDescriptor componentTypeDescriptor = parsingEngine.selectBestTypeDescriptor(componentType);

        for(Object id: multiplexedIds) {
            Object value = null;
            if(componentTypeDescriptor.createsTypeDefault(source))
                value = componentTypeDescriptor.createTypeDefault(parsingEngine,source,componentType);
            multiplexedMapping.put(id,value);
        }
        return multiplexedMapping;
    }

    @Override
    public String typeDefaultDocString(ArgumentSource source) {
        return "None";
    }

    @Override
    public Object parse(ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches) {
        if(multiplexedIds == null)
            throw new ReviewedGATKException("Cannot directly parse a MultiplexArgumentTypeDescriptor; must create a derivative type descriptor first.");

        Map<Object,Object> multiplexedMapping = new HashMap<Object,Object>();

        Class componentType = makeRawTypeIfNecessary(getCollectionComponentType(source.field));


        for(Object id: multiplexedIds) {
            Object value = parsingEngine.selectBestTypeDescriptor(componentType).parse(parsingEngine,source,componentType,matches.transform(multiplexer,id));
            multiplexedMapping.put(id,value);
        }

        parsingEngine.addTags(multiplexedMapping,getArgumentTags(matches));

        return multiplexedMapping;
    }

    public MultiplexArgumentTypeDescriptor createCustomTypeDescriptor(ParsingEngine parsingEngine,ArgumentSource dependentArgument,Object containingObject) {
        String[] sourceFields = dependentArgument.field.getAnnotation(Multiplex.class).arguments();

        List<ArgumentSource> allSources = parsingEngine.extractArgumentSources(containingObject.getClass());
        Class[] sourceTypes = new Class[sourceFields.length];
        Object[] sourceValues = new Object[sourceFields.length];
        int currentField = 0;

        for(String sourceField: sourceFields) {
            boolean fieldFound = false;
            for(ArgumentSource source: allSources) {
                if(!source.field.getName().equals(sourceField))
                    continue;
                if(source.field.isAnnotationPresent(Multiplex.class))
                    throw new ReviewedGATKException("Command-line arguments can only depend on independent fields");
                sourceTypes[currentField] = source.field.getType();
                sourceValues[currentField] = JVMUtils.getFieldValue(source.field,containingObject);
                currentField++;
                fieldFound = true;
            }
            if(!fieldFound)
                throw new ReviewedGATKException(String.format("Unable to find source field %s, referred to by dependent field %s",sourceField,dependentArgument.field.getName()));
        }

        Class<? extends Multiplexer> multiplexerType = dependentArgument.field.getAnnotation(Multiplex.class).value();
        Constructor<? extends Multiplexer> multiplexerConstructor;
        try {
            multiplexerConstructor = multiplexerType.getConstructor(sourceTypes);
            multiplexerConstructor.setAccessible(true);
        }
        catch(NoSuchMethodException ex) {
            throw new ReviewedGATKException(String.format("Unable to find constructor for class %s with parameters %s",multiplexerType.getName(),Arrays.deepToString(sourceFields)),ex);
        }

        Multiplexer multiplexer;
        try {
            multiplexer = multiplexerConstructor.newInstance(sourceValues);
        }
        catch(IllegalAccessException ex) {
            throw new ReviewedGATKException(String.format("Constructor for class %s with parameters %s is inaccessible",multiplexerType.getName(),Arrays.deepToString(sourceFields)),ex);
        }
        catch(InstantiationException ex) {
            throw new ReviewedGATKException(String.format("Can't create class %s with parameters %s",multiplexerType.getName(),Arrays.deepToString(sourceFields)),ex);
        }
        catch(InvocationTargetException ex) {
            throw new ReviewedGATKException(String.format("Can't invoke constructor of class %s with parameters %s",multiplexerType.getName(),Arrays.deepToString(sourceFields)),ex);
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
    protected Type getCollectionComponentType( Field field ) {
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
