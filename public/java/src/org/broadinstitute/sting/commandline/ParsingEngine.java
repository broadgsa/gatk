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

import com.google.java.contract.Requires;
import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.ApplicationDetails;
import org.broadinstitute.sting.utils.help.HelpFormatter;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.*;

/**
 * A parser for Sting command-line arguments.
 */
public class ParsingEngine {
    /**
     * The loaded argument sources along with their back definitions.
     */
    private Map<ArgumentDefinition,ArgumentSource> argumentSourcesByDefinition = new HashMap<ArgumentDefinition,ArgumentSource>();

    /**
     * A list of defined arguments against which command lines are matched.
     * Package protected for testing access.
     */
    public ArgumentDefinitions argumentDefinitions = new ArgumentDefinitions();

    /**
     * A list of matches from defined arguments to command-line text.
     * Indicates as best as possible where command-line text remains unmatched
     * to existing arguments.
     */
    ArgumentMatches argumentMatches = null;

    /**
     * Techniques for parsing and for argument lookup.
     */
    private List<ParsingMethod> parsingMethods = new ArrayList<ParsingMethod>();

    /**
     * All of the RodBinding objects we've seen while parsing
     */
    private List<RodBinding> rodBindings = new ArrayList<RodBinding>();

    /**
     * Class reference to the different types of descriptors that the create method can create.
     * The type of set used must be ordered (but not necessarily sorted).
     */
    private static final Set<ArgumentTypeDescriptor> STANDARD_ARGUMENT_TYPE_DESCRIPTORS = new LinkedHashSet<ArgumentTypeDescriptor>( Arrays.asList(new SimpleArgumentTypeDescriptor(),
            new IntervalBindingArgumentTypeDescriptor(),
            new RodBindingArgumentTypeDescriptor(),
            new CompoundArgumentTypeDescriptor(),
            new MultiplexArgumentTypeDescriptor()) );

    private Set<ArgumentTypeDescriptor> argumentTypeDescriptors = new LinkedHashSet<ArgumentTypeDescriptor>();

    /**
     * List of tags associated with the given instantiation of the command-line argument.
     */
    private final Map<Object,Tags> tags = new IdentityHashMap<Object,Tags>();    

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(ParsingEngine.class);

    public ParsingEngine( CommandLineProgram clp ) {
        RodBinding.resetNameCounter();
        parsingMethods.add( ParsingMethod.FullNameParsingMethod );
        parsingMethods.add( ParsingMethod.ShortNameParsingMethod );

        // Order matters here!  Make sure the clp's new type descriptors go in before the original type descriptors.
        if(clp != null)
            argumentTypeDescriptors.addAll(clp.getArgumentTypeDescriptors());
        argumentTypeDescriptors.addAll(STANDARD_ARGUMENT_TYPE_DESCRIPTORS);

        addArgumentSource(ParsingEngineArgumentFiles.class);
    }

    /**
     * Add a main argument source.  Argument sources are expected to have
     * any number of fields with an @Argument annotation attached.
     * @param source     An argument source from which to extract command-line arguments.
     */
    public void addArgumentSource( Class source ) {
        addArgumentSource(null, source);
    }

    /**
     * Add an argument source.  Argument sources are expected to have
     * any number of fields with an @Argument annotation attached.
     * @param sourceName name for this argument source.  'Null' indicates that this source should be treated
     *                   as the main module.
     * @param sourceClass A class containing argument sources from which to extract command-line arguments.
     */
    public void addArgumentSource( String sourceName, Class sourceClass ) {
        List<ArgumentDefinition> argumentsFromSource = new ArrayList<ArgumentDefinition>();
        for( ArgumentSource argumentSource: extractArgumentSources(sourceClass) ) {
            List<ArgumentDefinition> argumentDefinitions = argumentSource.createArgumentDefinitions();
            for(ArgumentDefinition argumentDefinition: argumentDefinitions) {
                argumentSourcesByDefinition.put(argumentDefinition,argumentSource);
                argumentsFromSource.add( argumentDefinition );
            }
        }
        argumentDefinitions.add( new ArgumentDefinitionGroup(sourceName, argumentsFromSource) );
    }

    /**
     * Do a cursory search to see if an argument with the given name is present.
     * @param argumentFullName full name of the argument.
     * @return True if the argument is present.  False otherwise.
     */
    public boolean isArgumentPresent( String argumentFullName ) {
        ArgumentDefinition definition =
                argumentDefinitions.findArgumentDefinition(argumentFullName,ArgumentDefinitions.FullNameDefinitionMatcher);
        return argumentMatches.hasMatch(definition);

    }

    /**
     * Parse the given set of command-line arguments, returning
     * an ArgumentMatches object describing the best fit of these
     * command-line arguments to the arguments that are actually
     * required.
     * @param tokens Tokens passed on the command line.
     * @return The parsed arguments by file.
     */
    public SortedMap<ArgumentMatchSource, List<String>> parse( String[] tokens ) {
        argumentMatches = new ArgumentMatches();
        SortedMap<ArgumentMatchSource, List<String>> parsedArgs = new TreeMap<ArgumentMatchSource, List<String>>();

        List<String> cmdLineTokens = Arrays.asList(tokens);
        parse(ArgumentMatchSource.COMMAND_LINE, cmdLineTokens, argumentMatches, parsedArgs);

        ParsingEngineArgumentFiles argumentFiles = new ParsingEngineArgumentFiles();

        // Load the arguments ONLY into the argument files.
        // Validation may optionally run on the rest of the arguments.
        loadArgumentsIntoObject(argumentFiles);

        for (File file: argumentFiles.files) {
            List<String> fileTokens = getArguments(file);
            parse(new ArgumentMatchSource(file), fileTokens, argumentMatches, parsedArgs);
        }

        return parsedArgs;
    }

    private void parse(ArgumentMatchSource matchSource, List<String> tokens,
                       ArgumentMatches argumentMatches, SortedMap<ArgumentMatchSource, List<String>> parsedArgs) {
        ArgumentMatchSite lastArgumentMatchSite = new ArgumentMatchSite(matchSource, -1);

        int i = 0;
        for (String token: tokens) {
            // If the token is of argument form, parse it into its own argument match.
            // Otherwise, pair it with the most recently used argument discovered.
            ArgumentMatchSite site = new ArgumentMatchSite(matchSource, i);
            if( isArgumentForm(token) ) {
                ArgumentMatch argumentMatch = parseArgument( token, site );
                if( argumentMatch != null ) {
                    argumentMatches.mergeInto( argumentMatch );
                    lastArgumentMatchSite = site;
                }
            }
            else {
                if( argumentMatches.hasMatch(lastArgumentMatchSite) &&
                    !argumentMatches.getMatch(lastArgumentMatchSite).hasValueAtSite(lastArgumentMatchSite))
                    argumentMatches.getMatch(lastArgumentMatchSite).addValue( lastArgumentMatchSite, token );
                else
                    argumentMatches.MissingArgument.addValue( site, token );

            }
            i++;
        }

        parsedArgs.put(matchSource, tokens);
    }

    private List<String> getArguments(File file) {
        try {
            if (file.getAbsolutePath().endsWith(".list")) {
                return getListArguments(file);
            }
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
        throw new UserException.CouldNotReadInputFile(file, "file extension is not .list");
    }

    private List<String> getListArguments(File file) throws IOException {
        ArrayList<String> argsList = new ArrayList<String>();
        for (String line: FileUtils.readLines(file))
            argsList.addAll(Arrays.asList(Utils.escapeExpressions(line)));
        return argsList;
    }

    public enum ValidationType { MissingRequiredArgument,
                                 InvalidArgument,
                                 InvalidArgumentValue,
                                 ValueMissingArgument,
                                 TooManyValuesForArgument,
                                 MutuallyExclusive }

    /**
     * Validates the list of command-line argument matches.
     */
    public void validate() {
        validate( EnumSet.noneOf(ValidationType.class) );
    }

    /**
     * Validates the list of command-line argument matches.  On failure throws an exception with detailed info about the
     * particular failures.  Takes an EnumSet indicating which validation checks to skip.
     * @param skipValidationOf List of validation checks to skip.
     */
    public void validate( EnumSet<ValidationType> skipValidationOf ) {
        // Find missing required arguments.
        if( !skipValidationOf.contains(ValidationType.MissingRequiredArgument) ) {
            Collection<ArgumentDefinition> requiredArguments =
                    argumentDefinitions.findArgumentDefinitions( true, ArgumentDefinitions.RequiredDefinitionMatcher );
            Collection<ArgumentDefinition> missingArguments = new ArrayList<ArgumentDefinition>();
            for( ArgumentDefinition requiredArgument: requiredArguments ) {
                if( !argumentMatches.hasMatch(requiredArgument) )
                    missingArguments.add( requiredArgument );
            }

            if( missingArguments.size() > 0 )
                throw new MissingArgumentException( missingArguments );
        }

        // Find invalid arguments.  Invalid arguments will have a null argument definition.
        if( !skipValidationOf.contains(ValidationType.InvalidArgument) ) {
            ArgumentMatches invalidArguments = argumentMatches.findUnmatched();
            if( invalidArguments.size() > 0 )
                throw new InvalidArgumentException( invalidArguments );
        }

        // Find invalid argument values -- invalid arguments are either completely missing or fail the specified 'validation' regular expression.
        if( !skipValidationOf.contains(ValidationType.InvalidArgumentValue) ) {
            Collection<ArgumentDefinition> verifiableArguments = 
                    argumentDefinitions.findArgumentDefinitions( null, ArgumentDefinitions.VerifiableDefinitionMatcher );
            Collection<Pair<ArgumentDefinition,String>> invalidValues = new ArrayList<Pair<ArgumentDefinition,String>>();
            for( ArgumentDefinition verifiableArgument: verifiableArguments ) {
                ArgumentMatches verifiableMatches = argumentMatches.findMatches( verifiableArgument );
                // Check to see whether an argument value was specified.  Argument values must be provided
                // when the argument name is specified and the argument is not a flag type.
                for(ArgumentMatch verifiableMatch: verifiableMatches) {
                    ArgumentSource argumentSource = argumentSourcesByDefinition.get(verifiableArgument);
                    if(verifiableMatch.values().size() == 0 && !verifiableArgument.isFlag && argumentSource.createsTypeDefault())
                        invalidValues.add(new Pair<ArgumentDefinition,String>(verifiableArgument,null));
                }

                // Ensure that the field contents meet the validation criteria specified by the regular expression.
                for( ArgumentMatch verifiableMatch: verifiableMatches ) {
                    for( String value: verifiableMatch.values() ) {
                        if( verifiableArgument.validation != null && !value.matches(verifiableArgument.validation) )
                            invalidValues.add( new Pair<ArgumentDefinition,String>(verifiableArgument, value) );
                    }
                }
            }

            if( invalidValues.size() > 0 )
                throw new InvalidArgumentValueException( invalidValues );
        }

        // Find values without an associated mate.
        if( !skipValidationOf.contains(ValidationType.ValueMissingArgument) ) {
            if( argumentMatches.MissingArgument.values().size() > 0 )
                throw new UnmatchedArgumentException( argumentMatches.MissingArgument );
        }

        // Find arguments with too many values.
        if( !skipValidationOf.contains(ValidationType.TooManyValuesForArgument)) {
            Collection<ArgumentMatch> overvaluedArguments = new ArrayList<ArgumentMatch>();
            for( ArgumentMatch argumentMatch: argumentMatches.findSuccessfulMatches() ) {
                // Warning: assumes that definition is not null (asserted by checks above).
                if( !argumentMatch.definition.isMultiValued && argumentMatch.values().size() > 1 )
                    overvaluedArguments.add(argumentMatch);
            }

            if( !overvaluedArguments.isEmpty() )
                throw new TooManyValuesForArgumentException(overvaluedArguments);
        }

        // Find sets of options that are supposed to be mutually exclusive.
        if( !skipValidationOf.contains(ValidationType.MutuallyExclusive)) {
            Collection<Pair<ArgumentMatch,ArgumentMatch>> invalidPairs = new ArrayList<Pair<ArgumentMatch,ArgumentMatch>>();
            for( ArgumentMatch argumentMatch: argumentMatches.findSuccessfulMatches() ) {
                if( argumentMatch.definition.exclusiveOf != null ) {
                    for( ArgumentMatch conflictingMatch: argumentMatches.findSuccessfulMatches() ) {
                        // Skip over the current element.
                        if( argumentMatch == conflictingMatch )
                            continue;
                        if( argumentMatch.definition.exclusiveOf.equals(conflictingMatch.definition.fullName) ||
                            argumentMatch.definition.exclusiveOf.equals(conflictingMatch.definition.shortName))
                            invalidPairs.add( new Pair<ArgumentMatch,ArgumentMatch>(argumentMatch, conflictingMatch) );
                    }
                }
            }

            if( !invalidPairs.isEmpty() )
                throw new ArgumentsAreMutuallyExclusiveException( invalidPairs );
        }
    }

    /**
     * Loads a set of matched command-line arguments into the given object.
     * @param object Object into which to add arguments.
     */
    public void loadArgumentsIntoObject( Object object ) {
        List<ArgumentSource> argumentSources = extractArgumentSources(object.getClass());

        List<ArgumentSource> dependentArguments = new ArrayList<ArgumentSource>();

        for( ArgumentSource argumentSource: argumentSources ) {
            if(argumentSource.isDeprecated() && argumentMatches.findMatches(this,argumentSource).size() > 0)
                notifyDeprecatedCommandLineArgument(argumentSource);

            // If this argument source depends on other command-line arguments, skip it and make a note to process it later.
            if(argumentSource.isDependent()) {
                dependentArguments.add(argumentSource);
                continue;
            }
            loadValueIntoObject( argumentSource, object, argumentMatches.findMatches(this,argumentSource) );
        }

        for(ArgumentSource dependentArgument: dependentArguments) {
            MultiplexArgumentTypeDescriptor dependentDescriptor = dependentArgument.createDependentTypeDescriptor(this,object);
            ArgumentSource dependentSource = dependentArgument.copyWithCustomTypeDescriptor(dependentDescriptor);
            loadValueIntoObject(dependentSource,object,argumentMatches.findMatches(this,dependentSource));
        }
    }

    /**
     * Notify the user that tags have been created.
     * @param key The key created.
     * @param tags List of tags, or empty list if no tags are present.
     */
    public void addTags(Object key, final Tags tags) {
        this.tags.put(key,tags);        
    }

    /**
     * Gets the tags associated with a given object.
     * @param key Key for which to find a tag.
     * @return List of tags associated with this key.
     */
    public Tags getTags(Object key)  {
        if(!tags.containsKey(key))
            return new Tags();
        return tags.get(key);
    }

    /**
     * Add a RodBinding type argument to this parser.  Called during parsing to allow
     * us to track all of the RodBindings discovered in the command line.
     * @param rodBinding the rodbinding to add.  Must not be added twice
     */
    @Requires("rodBinding != null")
    public void addRodBinding(final RodBinding rodBinding) {
        rodBindings.add(rodBinding);
    }

    /**
     * Notify the user that a deprecated command-line argument has been used.
     * @param argumentSource Deprecated argument source specified by user.
     */
    private void notifyDeprecatedCommandLineArgument(ArgumentSource argumentSource) {
        // Grab the first argument definition and report that one as the failure.  Theoretically, we should notify of all failures.
        List<ArgumentDefinition> definitions = argumentSource.createArgumentDefinitions();
        if(definitions.size() < 1)
            throw new ReviewedStingException("Internal error.  Argument source creates no definitions.");
        ArgumentDefinition definition = definitions.get(0);
        throw new UserException.DeprecatedArgument(definition.fullName,definition.doc);
    }

    /**
     * Loads a single argument into the object and that objects children.
     * @param argumentMatches Argument matches to load into the object.
     * @param source Argument source to load into the object.
     * @param instance Object into which to inject the value.  The target might be in a container within the instance.
     */
    private void loadValueIntoObject( ArgumentSource source, Object instance, ArgumentMatches argumentMatches ) {
        // Nothing to load
        if( argumentMatches.size() == 0 && ! source.createsTypeDefault() )
            return;

        // Target instance into which to inject the value.
        Collection<Object> targets = findTargets( source, instance );

        // Abort if no home is found for the object.
        if( targets.size() == 0 )
            throw new ReviewedStingException("Internal command-line parser error: unable to find a home for argument matches " + argumentMatches);

        for( Object target: targets ) {
            Object value = (argumentMatches.size() != 0) ? source.parse(this,argumentMatches) : source.createTypeDefault(this);

            JVMUtils.setFieldValue(source.field,target,value);
        }
    }

    public Collection<RodBinding> getRodBindings() {
        return Collections.unmodifiableCollection(rodBindings);
    }

    /**
     * Gets a collection of the container instances of the given type stored within the given target.
     * @param source Argument source.
     * @param instance Container.
     * @return A collection of containers matching the given argument source.
     */
    private Collection<Object> findTargets(ArgumentSource source, Object instance) {
        LinkedHashSet<Object> targets = new LinkedHashSet<Object>();
        for( Class clazz = instance.getClass(); clazz != null; clazz = clazz.getSuperclass() ) {
            for( Field field: clazz.getDeclaredFields() ) {
                if( field.equals(source.field) ) {
                    targets.add(instance);
                } else if( field.isAnnotationPresent(ArgumentCollection.class) ) {
                    targets.addAll(findTargets(source, JVMUtils.getFieldValue(field, instance)));
                }
            }
        }
        return targets;
    }

    /**
     * Prints out the help associated with these command-line argument definitions.
     * @param applicationDetails Details about the specific GATK-based application being run.
     */
    public void printHelp( ApplicationDetails applicationDetails ) {
        new HelpFormatter().printHelp(applicationDetails,argumentDefinitions);
    }

    /**
     * Extract all the argument sources from a given object.
     * @param sourceClass class to act as sources for other arguments.
     * @return A list of sources associated with this object and its aggregated objects.
     */
    public List<ArgumentSource> extractArgumentSources(Class sourceClass) {
        return extractArgumentSources(sourceClass, new Field[0]);
    }

    /**
     * Fetch the best command-line argument descriptor for the given class.
     * @param type Class for which to specify a descriptor.
     * @return descriptor for the given type.
     */
    public ArgumentTypeDescriptor selectBestTypeDescriptor(Class type) {
        return ArgumentTypeDescriptor.selectBest(argumentTypeDescriptors,type);
    }

    private List<ArgumentSource> extractArgumentSources(Class sourceClass, Field[] parentFields) {
        // now simply call into the truly general routine extract argument bindings but with a null
        // object so bindings aren't computed
        Map<ArgumentSource, Object> bindings = extractArgumentBindings(null, sourceClass, parentFields);
        return new ArrayList<ArgumentSource>(bindings.keySet());
    }

    public Map<ArgumentSource, Object> extractArgumentBindings(Object obj) {
        if ( obj == null ) throw new IllegalArgumentException("Incoming object cannot be null");
        return extractArgumentBindings(obj, obj.getClass(), new Field[0]);
    }

    /**
     * Extract all the argument sources from a given object, along with their bindings if obj != null .
     * @param obj the object corresponding to the sourceClass
     * @param sourceClass class to act as sources for other arguments.
     * @param parentFields Parent Fields
     * @return A map of sources associated with this object and its aggregated objects and bindings to their bindings values
     */
    private Map<ArgumentSource, Object> extractArgumentBindings(Object obj, Class sourceClass, Field[] parentFields) {
        Map<ArgumentSource, Object> bindings = new LinkedHashMap<ArgumentSource, Object>();

        while( sourceClass != null ) {
            Field[] fields = sourceClass.getDeclaredFields();
            for( Field field: fields ) {
                if( ArgumentTypeDescriptor.isArgumentAnnotationPresent(field) ) {
                    Object val = obj != null ? JVMUtils.getFieldValue(field, obj) : null;
                    bindings.put( new ArgumentSource(parentFields, field, selectBestTypeDescriptor(field.getType())), val );
                }
                if( field.isAnnotationPresent(ArgumentCollection.class) ) {
                    Object val = obj != null ? JVMUtils.getFieldValue(field, obj) : null;
                    Field[] newParentFields = Arrays.copyOf(parentFields, parentFields.length + 1);
                    newParentFields[parentFields.length] = field;
                    bindings.putAll( extractArgumentBindings(val, field.getType(), newParentFields) );
                }
            }

            sourceClass = sourceClass.getSuperclass();
        }

        return bindings;
    }

    /**
     * Determines whether a token looks like the name of an argument.
     * @param token Token to inspect.  Can be surrounded by whitespace.
     * @return True if token is of short name form.
     */
    private boolean isArgumentForm( String token ) {
        for( ParsingMethod parsingMethod: parsingMethods ) {
            if( parsingMethod.matches(token) )
                return true;
        }

        return false;
    }

    /**
     * Parse a short name into an ArgumentMatch.
     * @param token The token to parse.  The token should pass the isLongArgumentForm test.
     * @param position The position of the token in question.
     * @return ArgumentMatch associated with this token, or null if no match exists.
     */    
    private ArgumentMatch parseArgument( String token, ArgumentMatchSite position ) {
        if( !isArgumentForm(token) )
            throw new IllegalArgumentException( "Token is not recognizable as an argument: " + token );

        for( ParsingMethod parsingMethod: parsingMethods ) {
            if( parsingMethod.matches( token ) )
                return parsingMethod.match( argumentDefinitions, token, position );
        }

        // No parse results found.
        return null;
    }
}

/**
 * An exception indicating that some required arguments are missing.
 */
class MissingArgumentException extends ArgumentException {
    public MissingArgumentException( Collection<ArgumentDefinition> missingArguments ) {
        super( formatArguments(missingArguments) );
    }

    private static String formatArguments( Collection<ArgumentDefinition> missingArguments ) {
        StringBuilder sb = new StringBuilder();
        for( ArgumentDefinition missingArgument: missingArguments ) {
            if( missingArgument.shortName != null )
                sb.append( String.format("%nArgument with name '--%s' (-%s) is missing.", missingArgument.fullName, missingArgument.shortName) );
            else
                sb.append( String.format("%nArgument with name '--%s' is missing.", missingArgument.fullName) );
        }
        return sb.toString();
    }
}

/**
 * An exception for undefined arguments.
 */
class InvalidArgumentException extends ArgumentException {
    public InvalidArgumentException( ArgumentMatches invalidArguments ) {
        super( formatArguments(invalidArguments) );
    }

    private static String formatArguments( ArgumentMatches invalidArguments ) {
        StringBuilder sb = new StringBuilder();
        for( ArgumentMatch invalidArgument: invalidArguments )
            sb.append( String.format("%nArgument with name '%s' isn't defined.", invalidArgument.label) );
        return sb.toString();
    }
}

/**
 * An exception for values whose format is invalid.
 */
class InvalidArgumentValueException extends ArgumentException {
    public InvalidArgumentValueException( Collection<Pair<ArgumentDefinition,String>> invalidArgumentValues ) {
        super( formatArguments(invalidArgumentValues) );
    }

    private static String formatArguments( Collection<Pair<ArgumentDefinition,String>> invalidArgumentValues ) {
        StringBuilder sb = new StringBuilder();
        for( Pair<ArgumentDefinition,String> invalidValue: invalidArgumentValues ) {
            if(invalidValue.getSecond() == null)
                sb.append( String.format("%nArgument '--%s' requires a value but none was provided",
                                         invalidValue.first.fullName) );
            else
                sb.append( String.format("%nArgument '--%s' has value of incorrect format: %s (should match %s)",
                        invalidValue.first.fullName,
                        invalidValue.second,
                        invalidValue.first.validation) );
        }
        return sb.toString();
    }
}


/**
 * An exception for values that can't be mated with any argument.
 */
class UnmatchedArgumentException extends ArgumentException {
    public UnmatchedArgumentException( ArgumentMatch invalidValues ) {
        super( formatArguments(invalidValues) );
    }

    private static String formatArguments( ArgumentMatch invalidValues ) {
        StringBuilder sb = new StringBuilder();
        for( ArgumentMatchSite site: invalidValues.sites.keySet() )
            for( String value: invalidValues.sites.get(site) ) {
                switch (site.getSource().getType()) {
                    case CommandLine:
                        sb.append( String.format("%nInvalid argument value '%s' at position %d.",
                                value, site.getIndex()) );
                        break;
                    case File:
                        sb.append( String.format("%nInvalid argument value '%s' in file %s at position %d.",
                                value, site.getSource().getFile().getAbsolutePath(), site.getIndex()) );
                        break;
                    default:
                        throw new RuntimeException( String.format("Unexpected argument match source type: %s",
                                site.getSource().getType()));
                }
                if(value != null && Utils.dupString(' ',value.length()).equals(value))
                    sb.append("  Please make sure any line continuation backslashes on your command line are not followed by whitespace.");
            }
        return sb.toString();
    }
}

/**
 * An exception indicating that too many values have been provided for the given argument.
 */
class TooManyValuesForArgumentException extends ArgumentException {
    public TooManyValuesForArgumentException( Collection<ArgumentMatch> arguments ) {
        super( formatArguments(arguments) );
    }

    private static String formatArguments( Collection<ArgumentMatch> arguments ) {
        StringBuilder sb = new StringBuilder();
        for( ArgumentMatch argument: arguments )
            sb.append( String.format("%nArgument '%s' has too many values: %s.", argument.label, Arrays.deepToString(argument.values().toArray())) );
        return sb.toString();
    }
}

/**
 * An exception indicating that mutually exclusive options have been passed in the same command line.
 */
class ArgumentsAreMutuallyExclusiveException extends ArgumentException {
    public ArgumentsAreMutuallyExclusiveException( Collection<Pair<ArgumentMatch,ArgumentMatch>> arguments ) {
        super( formatArguments(arguments) );
    }

    private static String formatArguments( Collection<Pair<ArgumentMatch,ArgumentMatch>> arguments ) {
        StringBuilder sb = new StringBuilder();
        for( Pair<ArgumentMatch,ArgumentMatch> argument: arguments )
            sb.append( String.format("%nArguments '%s' and '%s' are mutually exclusive.", argument.first.definition.fullName, argument.second.definition.fullName ) );
        return sb.toString();
    }

}


/**
 * An exception for when an argument doesn't match an of the enumerated options for that var type
 */
class UnknownEnumeratedValueException extends ArgumentException {
    public UnknownEnumeratedValueException(ArgumentDefinition definition, String argumentPassed) {
        super( formatArguments(definition,argumentPassed) );
    }

    private static String formatArguments(ArgumentDefinition definition, String argumentPassed) {
        return String.format("Invalid value %s specified for argument %s; valid options are (%s).", argumentPassed, definition.fullName, Utils.join(",",definition.validOptions));
    }
}

/**
 * Container class to store the list of argument files.
 * The files will be parsed after the command line arguments.
 */
class ParsingEngineArgumentFiles {
    @Argument(fullName = "arg_file", shortName = "args", doc = "Reads arguments from the specified file", required = false)
    public List<File> files = new ArrayList<File>();
}
