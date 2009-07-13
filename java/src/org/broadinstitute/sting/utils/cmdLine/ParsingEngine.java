package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;
import org.apache.log4j.Logger;

import java.lang.reflect.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: May 3, 2009
 * Time: 4:35:25 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A parser for Sting command-line arguments.
 */
public class ParsingEngine {
    /**
     * A collection of all the source fields which define command-line arguments.
     */
    List<ArgumentSource> argumentSources = new ArrayList<ArgumentSource>();

    /**
     * A list of defined arguments against which command lines are matched.
     * Package protected for testing access.
     */
    ArgumentDefinitions argumentDefinitions = new ArgumentDefinitions();

    /**
     * A list of matches from defined arguments to command-line text.
     * Indicates as best as possible where command-line text remains unmatched
     * to existing arguments.
     */
    ArgumentMatches argumentMatches = null;

    /**
     * Stores a custom argument factory for building out arguments of which only
     * subclasses of CommandLineProgram should be aware.
     */
    ArgumentFactory customArgumentFactory = null;


    /**
     * Techniques for parsing and for argument lookup.
     */
    private List<ParsingMethod> parsingMethods = new ArrayList<ParsingMethod>();

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(ParsingEngine.class);

    public ParsingEngine( ArgumentFactory customArgumentFactory ) {
        this.customArgumentFactory = customArgumentFactory;
        parsingMethods.add( ParsingMethod.FullNameParsingMethod );
        parsingMethods.add( ParsingMethod.ShortNameParsingMethod );
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
        for( ArgumentSource argumentSource: extractArgumentSources(sourceClass,true) )
            argumentsFromSource.add( new ArgumentDefinition(argumentSource) );
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
     * @return A object indicating which matches are best.  Might return
     *         an empty object, but will never return null.
     */
    public void parse( String[] tokens ) {
        argumentMatches = new ArgumentMatches();

        int lastArgumentMatchSite = -1;

        for( int i = 0; i < tokens.length; i++ ) {
            String token = tokens[i];
            // If the token is of argument form, parse it into its own argument match.
            // Otherwise, pair it with the most recently used argument discovered.
            if( isArgumentForm(token) ) {
                ArgumentMatch argumentMatch = parseArgument( token, i );
                if( argumentMatch != null ) {
                    argumentMatches.mergeInto( argumentMatch );
                    lastArgumentMatchSite = i;
                }
            }
            else {
                if( argumentMatches.hasMatch(lastArgumentMatchSite) &&
                    !argumentMatches.getMatch(lastArgumentMatchSite).hasValueAtSite(lastArgumentMatchSite))
                    argumentMatches.getMatch(lastArgumentMatchSite).addValue( lastArgumentMatchSite, token );
                else
                    argumentMatches.MissingArgument.addValue( i, token );

            }
        }
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
            Collection<ArgumentMatch> invalidArguments = argumentMatches.findUnmatched();
            if( invalidArguments.size() > 0 )
                throw new InvalidArgumentException( invalidArguments );
        }

        // Find invalid argument values (arguments that fail the regexp test.
        if( !skipValidationOf.contains(ValidationType.InvalidArgumentValue) ) {
            Collection<ArgumentDefinition> verifiableArguments = 
                    argumentDefinitions.findArgumentDefinitions( null, ArgumentDefinitions.VerifiableDefinitionMatcher );
            Collection<Pair<ArgumentDefinition,String>> invalidValues = new ArrayList<Pair<ArgumentDefinition,String>>();
            for( ArgumentDefinition verifiableArgument: verifiableArguments ) {
                Collection<ArgumentMatch> verifiableMatches = argumentMatches.findMatches( verifiableArgument );
                for( ArgumentMatch verifiableMatch: verifiableMatches ) {
                    for( String value: verifiableMatch.values() ) {
                        if( !value.matches(verifiableArgument.validation) )
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
                if( !argumentMatch.definition.source.isMultiValued() && argumentMatch.values().size() > 1 )
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
        // Get a list of argument sources, not including the children of this argument.  For now, skip loading
        // arguments into the object recursively.
        List<ArgumentSource> argumentSources = extractArgumentSources( object.getClass(), false );
        for( ArgumentSource argumentSource: argumentSources ) {
            Collection<ArgumentMatch> argumentsMatchingSource = argumentMatches.findMatches( argumentSource );
            if( argumentsMatchingSource.size() != 0 )
                loadMatchesIntoObject( argumentsMatchingSource, object );
        }
    }

    /**
     * Loads a single argument into the object.
     * @param argumentMatches Argument matches to load into the object.
     * @param object Target for the argument.
     */
    private void loadMatchesIntoObject( Collection<ArgumentMatch> argumentMatches, Object object ) {
        if( argumentMatches.size() > 1 )
            throw new StingException("Too many matches");

        ArgumentMatch match = argumentMatches.iterator().next();
        ArgumentDefinition definition = match.definition;

        // A null definition might be in the list if some invalid arguments were passed in but we
        // want to load in a subset of data for better error reporting.  Ignore null definitions.
        if( definition == null )
            return;

        if( definition.source.clazz.isAssignableFrom(object.getClass()) ) {
            if( !definition.source.isFlag() ) {
                String[] tokens = match.values().toArray(new String[0]);
                FieldParser fieldParser = FieldParser.create(definition.source.field);
                definition.source.setValue( object, fieldParser.parse(tokens) );
            }
            else
                definition.source.setValue( object, true );
        }
    }

    /**
     * Prints out the help associated with these command-line argument definitions.
     */
    public void printHelp( String runningInstructions ) {
        new HelpFormatter().printHelp(runningInstructions,argumentDefinitions);
    }

    /**
     * Extract all the argument sources from a given object.
     * @param sourceClass class to act as sources for other arguments.
     * @param recursive Whether to recursively look for argument collections and add their contents.
     * @return A list of sources associated with this object and its aggregated objects.
     */
    private List<ArgumentSource> extractArgumentSources( Class sourceClass, boolean recursive ) {
        List<ArgumentSource> argumentSources = new ArrayList<ArgumentSource>();
        while( sourceClass != null ) {
            Field[] fields = sourceClass.getDeclaredFields();
            for( Field field: fields ) {
                if( field.isAnnotationPresent(Argument.class) )
                    argumentSources.add( new ArgumentSource(sourceClass,field) );
                if( field.isAnnotationPresent(ArgumentCollection.class) && recursive )
                    argumentSources.addAll( extractArgumentSources(field.getType(),recursive) );
            }
            sourceClass = sourceClass.getSuperclass();
        }
        return argumentSources;
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
     * @return ArgumentMatch associated with this token, or null if no match exists.
     */    
    private ArgumentMatch parseArgument( String token, int position ) {
        if( !isArgumentForm(token) )
            throw new IllegalArgumentException( "Token is not recognizable as an argument: " + token );

        for( ParsingMethod parsingMethod: parsingMethods ) {
            if( parsingMethod.matches( token ) )
                return parsingMethod.match( argumentDefinitions, token, position );
        }

        // No parse results found.
        return null;
    }

    /**
     * Constructs a command-line argument given a string and field.
     * @param f Field type from which to infer the type.
     * @param strs Collection of parameter strings to parse.
     * @return Parsed object of the inferred type.
     */
    private Object constructFromString(Field f, List<String> strs) {
        FieldParser fieldParser = FieldParser.create(f);
        return fieldParser.parse( strs.toArray(new String[0]) );
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
    public InvalidArgumentException( Collection<ArgumentMatch> invalidArguments ) {
        super( formatArguments(invalidArguments) );
    }

    private static String formatArguments( Collection<ArgumentMatch> invalidArguments ) {
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
        for( int index: invalidValues.indices.keySet() )
            for( String value: invalidValues.indices.get(index) )
                sb.append( String.format("%nInvalid argument value '%s' at position %d.", value, index) );
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
            sb.append( String.format("%nArgument '%s' has to many values: %s.", argument.label, Arrays.deepToString(argument.values().toArray())) );
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