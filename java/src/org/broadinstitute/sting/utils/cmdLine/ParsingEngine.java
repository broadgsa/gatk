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
     */
    private ArgumentDefinitions argumentDefinitions = new ArgumentDefinitions();

    /**
     * Techniques for parsing and for argument lookup.
     */
    private List<ParsingMethod> parsingMethods = new ArrayList<ParsingMethod>();

    /**
     * our log, which we want to capture anything from org.broadinstitute.sting
     */
    protected static Logger logger = Logger.getLogger(ArgumentParser.class);    

    public ParsingEngine() {
        parsingMethods.add( new ParsingMethod(Pattern.compile("\\s*--([\\w\\.]+)\\s*"), new FullNameDefinitionMatcher()) );
        parsingMethods.add( new ParsingMethod(Pattern.compile("\\s*-([\\w\\.]+)\\s*"), new ShortNameDefinitionMatcher()) );
        parsingMethods.add( new ParsingMethod(Pattern.compile("\\s*-([\\w\\.])([\\w\\.]+)\\s*"), new ShortNameDefinitionMatcher()) );
    }

    /**
     * Add an argument source.  Argument sources are expected to have
     * any number of fields with an @Argument annotation attached.
     * @param sources A list of argument sources from which to extract
     *                command-line arguments.
     */
    public void addArgumentSources( Class... sources ) {
        for( Class source: sources ) {
            Field[] fields = source.getFields();
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

    /**
     * Validates the list of command-line argument matches.  On
     * failure ...TBD...
     */
    public void validate( ArgumentMatches argumentMatches ) {
        
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
                    if( !isArgumentBoolean(definition) )
                        definition.sourceField.set( object, constructFromString( definition.sourceField, match.values ) );
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

    private boolean isArgumentBoolean( ArgumentDefinition definition ) {
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
            if( parsingMethod.hasMatch( argumentDefinitions, token ) )
                return parsingMethod.findMatch( argumentDefinitions, token, position );
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
        ArgumentMatch lastMatched = null;

        for( int i = 0; i < tokens.length; i++ ) {
            if( argumentMatches.hasMatch(i) ) {
                lastMatched = argumentMatches.getMatch(i);
                continue;
            }
            
            lastMatched.addValue( tokens[i] );
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

        public ParsingMethod( Pattern pattern, DefinitionMatcher definitionMatcher ) {
            this.pattern = pattern;
            this.definitionMatcher = definitionMatcher;
        }

        public boolean hasMatch( ArgumentDefinitions definitions, String token ) {
            Matcher matcher = pattern.matcher(token);
            return matcher.matches() && definitionMatcher.get( definitions, matcher.group(1) ) != null;
        }

        public ArgumentMatch findMatch( ArgumentDefinitions definitions, String token, int position ) {
            Matcher matcher = pattern.matcher(token);

            // Didn't match?  Must be bad input.
            if( !matcher.matches() )
                throw new IllegalArgumentException( String.format("Unable to parse token %s with pattern %s", token, pattern.pattern()) );

            // If the argument is valid, parse out the argument and value (if present).
            String argument = matcher.group(1);
            String value = matcher.groupCount() > 1 ? matcher.group(2) : null;

            // Try to find a matching argument.  If found, label that as the match.
            ArgumentDefinition argumentDefinition = definitionMatcher.get( definitions, argument );

            if( argumentDefinition != null ) {
                ArgumentMatch argumentMatch = new ArgumentMatch( argumentDefinition, position );
                if( value != null )
                    argumentMatch.addValue( value );
                return argumentMatch;
            }

            throw new IllegalArgumentException( String.format("Unable to find match for token %s", token) );
        }
    }
}
