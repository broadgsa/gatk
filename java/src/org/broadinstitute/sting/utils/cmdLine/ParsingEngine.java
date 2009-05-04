package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;

import java.lang.reflect.Field;
import java.util.Set;

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
     * @param arguments Command-line arguments.
     * @return A object indicating which matches are best.  Might return
     *         an empty object, but will never return null.
     */
    public ArgumentMatches parse( String[] arguments ) {
        ArgumentMatches argumentMatches = new ArgumentMatches();

        for( int i = 0; i < arguments.length; i++ ) {
            String argument = arguments[i].trim();
            if( argument.startsWith("-") ) {
                String shortName = argument.substring(1);
                if( argumentDefinitions.hasArgumentWithShortName(shortName) ) {
                    ArgumentDefinition definition = argumentDefinitions.getArgumentWithShortName(shortName);
                    argumentMatches.add( definition, arguments[i+1].trim() );
                }
            }
        }

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
     * @param argumentMatches List of matches.
     */
    public void loadArgumentsIntoObject( Object object, ArgumentMatches matches ) {
        for( ArgumentMatch match: matches ) {
            if( object.getClass().equals(match.definition.sourceClass) ) {
                try {
                    match.definition.sourceField.set( object, match.value );
                }
                catch( IllegalAccessException ex ) {
                    //logger.fatal("processArgs: cannot convert field " + field.toString());
                    throw new StingException("processArgs: Failed conversion " + ex.getMessage(), ex);                    
                }
            }
        }
    }
}
