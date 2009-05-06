package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;
import org.apache.log4j.Logger;

import java.lang.reflect.Field;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Modifier;
import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
import java.util.Arrays;
import java.util.EnumSet;

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
     * A list of defined arguments against which command lines are matched.
     * Package protected for testing access.
     */
    ArgumentDefinitions argumentDefinitions = new ArgumentDefinitions();

    /**
     * Techniques for parsing and for argument lookup.
     */
    private List<ParsingMethod> parsingMethods = new ArrayList<ParsingMethod>();

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(ArgumentParser.class);

    public ParsingEngine() {
        parsingMethods.add( new ParsingMethod(Pattern.compile("\\s*--([\\w\\.]+)\\s*"), ArgumentDefinitions.FullNameDefinitionMatcher) );
        parsingMethods.add( new ParsingMethod(Pattern.compile("\\s*-([\\w\\.])([\\w\\.]*)\\s*"),
                            ArgumentDefinitions.ShortNameDefinitionMatcher,
                            ArgumentDefinitions.ShortNameAliasProvider) );
    }

    /**
     * Add an argument source.  Argument sources are expected to have
     * any number of fields with an @Argument annotation attached.
     * @param sources A list of argument sources from which to extract
     *                command-line arguments.
     */
    public void addArgumentSources( Class... sources ) {
        for( Class source: sources ) {
            Field[] fields = source.getDeclaredFields();
            for( Field field: fields ) {
                Argument argument = field.getAnnotation(Argument.class);                
                if(argument != null)
                    argumentDefinitions.add( argument, source, field );
            }
        }
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
    public ArgumentMatches parse( String[] tokens ) {
        ArgumentMatches argumentMatches = parseArguments( tokens );
        fitValuesToArguments( argumentMatches, tokens );
        return argumentMatches;
    }

    public enum ValidationType { MissingRequiredArgument,
                                 InvalidArgument,
                                 ValueMissingArgument,
                                 TooManyValuesForArgument };

    /**
     * Validates the list of command-line argument matches.
     * @param argumentMatches Matches to validate.
     */
    public void validate( ArgumentMatches argumentMatches ) {
        validate( argumentMatches, EnumSet.noneOf(ValidationType.class) );
    }

    /**
     * Validates the list of command-line argument matches.  On failure throws an exception with detailed info about the
     * particular failures.  Takes an EnumSet indicating which validation checks to skip.
     * @param argumentMatches Matches to validate.
     * @param skipValidationOf List of validation checks to skip.
     */
    public void validate( ArgumentMatches argumentMatches, EnumSet<ValidationType> skipValidationOf ) {
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
            Collection<ArgumentMatch> invalidArguments = argumentMatches.findMatches(null);
            if( invalidArguments.size() > 0 )
                throw new InvalidArgumentException( invalidArguments );
        }

        // Find values without an associated mate.
        if( !skipValidationOf.contains(ValidationType.ValueMissingArgument) ) {
            if( argumentMatches.MissingArgument.values().size() > 0 )
                throw new InvalidArgumentValueException( argumentMatches.MissingArgument );
        }

        // Find arguments with too many values.
        if( !skipValidationOf.contains(ValidationType.TooManyValuesForArgument)) {
            Collection<ArgumentMatch> overvaluedArguments = new ArrayList<ArgumentMatch>();
            for( ArgumentMatch argumentMatch: argumentMatches ) {
                // Warning: assumes that definition is not null (asserted by checks above).
                if( argumentMatch.definition != null &&
                        !argumentMatch.definition.isMultiValued() &&
                        argumentMatch.values().size() > 1 )
                    overvaluedArguments.add(argumentMatch);
            }

            if( !overvaluedArguments.isEmpty() )
                throw new TooManyValuesForArgumentException(overvaluedArguments);
        }
    }

    /**
     * Loads a set of matched command-line arguments into the given object.
     * @param object Object into which to add arguments.
     * @param matches List of matches.
     */
    public void loadArgumentsIntoObject( Object object, ArgumentMatches matches ) {
        for( ArgumentMatch match: matches ) {
            ArgumentDefinition definition = match.definition;
            if( object.getClass().equals(definition.sourceClass) ) {
                try {
                    if( !isArgumentFlag(definition) )
                        definition.sourceField.set( object, constructFromString( definition.sourceField, match.values() ) );
                    else
                        definition.sourceField.set( object, true );
                }
                catch( IllegalAccessException ex ) {
                    //logger.fatal("processArgs: cannot convert field " + field.toString());
                    throw new StingException("processArgs: Failed conversion " + ex.getMessage(), ex);                    
                }
            }
        }
    }

    /**
     * Returns true if the argument is a flag (a 0-valued argument).
     * @param definition Argument definition.
     * @return True if argument is a flag; false otherwise.
     */
    private boolean isArgumentFlag( ArgumentDefinition definition ) {
        return (definition.sourceField.getType() == Boolean.class) || (definition.sourceField.getType() == Boolean.TYPE);
    }

    /**
     * Determines whether a token looks like the name of an argument.
     * @param token Token to inspect.  Can be surrounded by whitespace.
     * @return True if token is of short name form.
     */
    private boolean isArgumentForm( String token ) {
        for( ParsingMethod parsingMethod: parsingMethods ) {
            if( parsingMethod.pattern.matcher(token).matches() )
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
            if( parsingMethod.matches( argumentDefinitions, token ) )
                return parsingMethod.match( argumentDefinitions, token, position );
        }

        // No parse results found.
        return null;
    }

    /**
     * Extracts the argument portions of the string and assemble them into a data structure.
     * @param tokens List of tokens from which to find arguments.
     * @return Set of argument matches.
     */
    private ArgumentMatches parseArguments( String[] tokens ) {
        ArgumentMatches argumentMatches = new ArgumentMatches();

        for( int i = 0; i < tokens.length; i++ ) {
            String token = tokens[i];
            if( isArgumentForm(token) ) {
                ArgumentMatch argumentMatch = parseArgument( token, i );
                if( argumentMatch != null )
                    argumentMatches.mergeInto( argumentMatch );
            }
        }

        return argumentMatches;
    }

    /**
     * Fit the options presented on the command line to the given arguments.
     * @param argumentMatches List of arguments already matched to data.
     * @param tokens The command-line input.
     */
    private void fitValuesToArguments( ArgumentMatches argumentMatches, String[] tokens ) {
        for( int i = 0; i < tokens.length; i++ ) {
            // If this is the site of a successfully matched argument, pass it over.
            if( argumentMatches.hasMatch(i) )
                continue;

            // tokens[i] must be an argument value.  Match it with the previous argument.
            String value = tokens[i];
            int argumentSite = i - 1;

            // If the argument is present and doesn't already have a value associated with the given site, add the value.
            if( argumentMatches.hasMatch(argumentSite) && !argumentMatches.getMatch(argumentSite).hasValueAtSite(argumentSite))
                argumentMatches.getMatch(argumentSite).addValue( argumentSite, value );
            else
                argumentMatches.MissingArgument.addValue( i, value );
        }
    }

    /**
     * Constructs a command-line argument given a string and field.
     * @param f Field type from which to infer the type.
     * @param strs Collection of parameter strings to parse.
     * @return Parsed object of the inferred type.
     */
    private Object constructFromString(Field f, List<String> strs) {
        Class type = f.getType();

        if( Collection.class.isAssignableFrom(type) ) {
            Collection collection = null;
            Class containedType = null;

            // If this is a parameterized collection, find the contained type.  If blow up if only one type exists.
            if( f.getGenericType() instanceof ParameterizedType) {
                ParameterizedType parameterizedType = (ParameterizedType)f.getGenericType();
                if( parameterizedType.getActualTypeArguments().length > 1 )
                    throw new IllegalArgumentException("Unable to determine collection type of field: " + f.toString());
                containedType = (Class)parameterizedType.getActualTypeArguments()[0];
            }
            else
                containedType = String.class;

            // If this is a generic interface, pick a concrete implementation to create and pass back.
            // Because of type erasure, don't worry about creating one of exactly the correct type.
            if( Modifier.isInterface(type.getModifiers()) || Modifier.isAbstract(type.getModifiers()) )
            {
                if( java.util.List.class.isAssignableFrom(type) ) type = ArrayList.class;
                else if( java.util.Queue.class.isAssignableFrom(type) ) type = java.util.ArrayDeque.class;
                else if( java.util.Set.class.isAssignableFrom(type) ) type = java.util.TreeSet.class;
            }

            try
            {
                collection = (Collection)type.newInstance();
            }
            catch( Exception ex ) {
                // Runtime exceptions are definitely unexpected parsing simple collection classes.
                throw new IllegalArgumentException(ex);
            }

            for( String str: strs )
                collection.add( constructSingleElement(f,containedType,str) );

            return collection;
        }
        else if( type.isArray() ) {
            Class containedType = type.getComponentType();

            Object arr = Array.newInstance(containedType,strs.size());
            for( int i = 0; i < strs.size(); i++ )
                Array.set( arr,i,constructSingleElement(f,containedType,strs.get(i)) );
            return arr;
        }
        else  {
            if( strs.size() != 1 )
                throw new IllegalArgumentException("Passed multiple arguments to an object expecting a single value.");
            return constructSingleElement(f,type,strs.get(0));
        }
    }

    /**
     * Builds a single element of the given type.
     * @param f Implies type of data to construct.
     * @param str String representation of data.
     * @return parsed form of String.
     */
    private Object constructSingleElement(Field f, Class type, String str) {
        // lets go through the types we support
        if (type == Boolean.TYPE) {
            boolean b = false;
            if (str.toLowerCase().equals("true")) {
                b = true;
            }
            Boolean bool = new Boolean(b);
            return bool;
        } else if (type == Integer.TYPE) {
            Integer in = Integer.valueOf(str);
            return in;
        } else if (type == Float.TYPE) {
            Float fl = Float.valueOf(str);
            return fl;
        }
        else {
            Constructor ctor = null;
            try {
                ctor = type.getConstructor(String.class);
                return ctor.newInstance(str);
            } catch (NoSuchMethodException e) {
                logger.fatal("constructFromString:NoSuchMethodException: cannot convert field " + f.toString());
                throw new RuntimeException("constructFromString:NoSuchMethodException: Failed conversion " + e.getMessage());
            } catch (IllegalAccessException e) {
                logger.fatal("constructFromString:IllegalAccessException: cannot convert field " + f.toString());
                throw new RuntimeException("constructFromString:IllegalAccessException: Failed conversion " + e.getMessage());
            } catch (InvocationTargetException e) {
                logger.fatal("constructFromString:InvocationTargetException: cannot convert field " + f.toString());
                throw new RuntimeException("constructFromString:InvocationTargetException: Failed conversion " + e.getMessage());
            } catch (InstantiationException e) {
                logger.fatal("constructFromString:InstantiationException: cannot convert field " + f.toString());
                throw new RuntimeException("constructFromString:InstantiationException: Failed conversion " + e.getMessage());
            }
        }
    }


    
    /**
     * Holds a pattern, along with how to get to the argument definitions that could match that pattern.
     */
    private class ParsingMethod {
        public final Pattern pattern;
        public final DefinitionMatcher definitionMatcher;
        public final AliasProvider aliasProvider;

        public ParsingMethod( Pattern pattern, DefinitionMatcher definitionMatcher ) {
            this( pattern, definitionMatcher, null );
        }

        public ParsingMethod( Pattern pattern, DefinitionMatcher definitionMatcher, AliasProvider aliasProvider ) {
            this.pattern = pattern;
            this.definitionMatcher = definitionMatcher;
            this.aliasProvider = aliasProvider;
        }

        public boolean matches( ArgumentDefinitions definitions, String token ) {
            Matcher matcher = pattern.matcher(token);
            return matcher.matches();
        }

        public ArgumentMatch match( ArgumentDefinitions definitions, String token, int position ) {
            Matcher matcher = pattern.matcher(token);

            // Didn't match?  Must be bad input.
            if( !matcher.matches() )
                throw new IllegalArgumentException( String.format("Unable to parse token %s with pattern %s", token, pattern.pattern()) );

            // If the argument is valid, parse out the argument and value (if present).
            String argument = matcher.group(1);
            String value = null;
            if( matcher.groupCount() > 1 && matcher.group(2).trim().length() > 0)
                value = matcher.group(2).trim();

            // If an alias provider has been provided, determine the possible list of argument names that this
            // argument / value pair can represent.
            ArgumentDefinition bestMatchArgumentDefinition = null;
            if( aliasProvider != null ) {
                List<String> aliases = aliasProvider.getAliases( argument, value );
                String bestAlias = null;

                for( String alias: aliases ) {
                    if( definitions.findArgumentDefinition(alias,definitionMatcher) != null ) {
                        bestAlias = alias;
                        bestMatchArgumentDefinition = definitions.findArgumentDefinition(alias,definitionMatcher);
                        break;
                    }
                }

                // Couldn't find anything appropriate?  The aliases should be in best-to-worst order, so
                if( bestAlias == null ) {
                    bestAlias = aliases.get(0);
                }

                if( aliasProvider.doesAliasConsumeValue(bestAlias,argument,value) ) value = null;
                argument = bestAlias;                
            }
            else
                bestMatchArgumentDefinition = definitions.findArgumentDefinition( argument, definitionMatcher );                

            // Try to find a matching argument.  If found, label that as the match.  If not found, add the argument
                // with a null definition.
            ArgumentMatch argumentMatch = new ArgumentMatch( argument, bestMatchArgumentDefinition, position );
            if( value != null )
                argumentMatch.addValue( position, value );
            return argumentMatch;
        }
    }
}

/**
 * An exception indicating that some required arguments are missing.
 */
class MissingArgumentException extends StingException {
    public MissingArgumentException( Collection<ArgumentDefinition> missingArguments ) {
        super( formatArguments(missingArguments) );
    }

    private static String formatArguments( Collection<ArgumentDefinition> missingArguments ) {
        StringBuilder sb = new StringBuilder();
        for( ArgumentDefinition missingArgument: missingArguments )
            sb.append( String.format("Argument with name '%s' is missing.", missingArgument.fullName) );
        return sb.toString();
    }
}

/**
 * An exception for undefined arguments.
 */
class InvalidArgumentException extends StingException {
    public InvalidArgumentException( Collection<ArgumentMatch> invalidArguments ) {
        super( formatArguments(invalidArguments) );
    }

    private static String formatArguments( Collection<ArgumentMatch> invalidArguments ) {
        StringBuilder sb = new StringBuilder();
        for( ArgumentMatch invalidArgument: invalidArguments )
            sb.append( String.format("Argument with name '%s' isn't defined.", invalidArgument.label) );
        return sb.toString();
    }
}

/**
 * An exception for values that can't be mated with any argument.
 */
class InvalidArgumentValueException extends StingException {
    public InvalidArgumentValueException( ArgumentMatch invalidValues ) {
        super( formatArguments(invalidValues) );
    }

    private static String formatArguments( ArgumentMatch invalidValues ) {
        StringBuilder sb = new StringBuilder();
        for( int index: invalidValues.indices.keySet() )
            for( String value: invalidValues.indices.get(index) )
                sb.append( String.format("Invalid argument value '%s' at position %d", value, index) );
        return sb.toString();
    }
}

/**
 * An exception indicating that too many values have been provided for the given argument.
 */
class TooManyValuesForArgumentException extends StingException {
    public TooManyValuesForArgumentException( Collection<ArgumentMatch> arguments ) {
        super( formatArguments(arguments) );
    }

    private static String formatArguments( Collection<ArgumentMatch> arguments ) {
        StringBuilder sb = new StringBuilder();
        for( ArgumentMatch argument: arguments )
            sb.append( String.format("Argument with name '%s' has to many values: %s", argument.label, Arrays.deepToString(argument.values().toArray())) );
        return sb.toString();
    }
}