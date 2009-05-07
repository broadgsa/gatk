package org.broadinstitute.sting.utils.cmdLine;

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.List;
import java.util.ArrayList;

/**
 * Holds a pattern, along with how to get to the argument definitions that could match that pattern.
 */
interface ParsingMethod {
    /**
     * Can the given token be parsed by this parsing method?
     * @param token Token to validate.
     * @return True if the given token matches.
     */
    public abstract boolean matches( String token );

    /**
     * Find the best match for a given token at a given position from among the provided
     * argument definitions.
     * @param definitions List of argument definitions.
     * @param token The token from the command line to match.  Should be validated using
     *              ParsingMethod's matches() tester.
     * @param position Position at which this command-line argument occurs.  Will be used
     *                 for validation later.
     * @return An argument match.  Definition field will be populated if a match was found or
     *         empty if no appropriate definition could be found. 
     */
    public abstract ArgumentMatch match( ArgumentDefinitions definitions, String token, int position );
}

/**
 * Instructions for how to parse a command-line argument passed by full name into a match.
 */
class FullNameParsingMethod implements ParsingMethod {
    private static final Pattern pattern = Pattern.compile("\\s*--([\\w\\.\\-]+)\\s*");
    private static final DefinitionMatcher definitionMatcher = ArgumentDefinitions.FullNameDefinitionMatcher;

    public boolean matches( String token ) {
        Matcher matcher = pattern.matcher(token);
        return matcher.matches();
    }

    public ArgumentMatch match( ArgumentDefinitions definitions, String token, int position ) {
        // If the argument is valid, parse out the argument.
        Matcher matcher = pattern.matcher(token);

        // Didn't match?  Must be bad input.
        if( !matcher.matches() )
            throw new IllegalArgumentException( String.format("Unable to parse token %s with pattern %s", token, pattern.pattern()) );

        String argument = matcher.group(1).trim();

        // Find the most appropriate argument definition for the given argument.
        ArgumentDefinition argumentDefinition = definitions.findArgumentDefinition( argument, definitionMatcher );

        // Try to find a matching argument.  If found, label that as the match.  If not found, add the argument
        // with a null definition.
        ArgumentMatch argumentMatch = new ArgumentMatch( argument, argumentDefinition, position );

        return argumentMatch;
    }
}

/**
 * Instructions for how to parse a command-line argument passed by short name into a match.
 */
class ShortNameParsingMethod implements ParsingMethod {
    private static final Pattern standalonePattern = Pattern.compile("\\s*-([\\w\\-]+)\\s*");
    private static final Pattern embeddedValuePattern = Pattern.compile("\\s*-([\\w\\.])([\\w/:\\.\\-]+)\\s*");
    private static final DefinitionMatcher definitionMatcher = ArgumentDefinitions.ShortNameDefinitionMatcher;

    public boolean matches( String token ) {
        return standalonePattern.matcher(token).matches() || embeddedValuePattern.matcher(token).matches();
    }    

    public ArgumentMatch match( ArgumentDefinitions definitions, String token, int position ) {
        // Didn't match?  Must be bad input.
        if( !matches(token) )
            throw new IllegalArgumentException( String.format("Unable to parse token %s with pattern %s", token, embeddedValuePattern.pattern()) );

        // Build the best possible standalone match given the available data.
        Matcher standaloneMatcher = standalonePattern.matcher(token);
        ArgumentMatch standaloneMatch = null;

        if( standaloneMatcher.matches() ) {
            String argument = standaloneMatcher.group(1).trim();

            ArgumentDefinition argumentDefinition = definitions.findArgumentDefinition(argument,definitionMatcher);
            standaloneMatch = new ArgumentMatch( argument, argumentDefinition, position );
        }

        // Build the best possible embedded value match given the available data.
        Matcher embeddedValueMatcher = embeddedValuePattern.matcher(token);
        ArgumentMatch embeddedValueMatch = null;

        if( embeddedValueMatcher.matches() ) {
            String argument = embeddedValueMatcher.group(1).trim();
            String value = embeddedValueMatcher.group(2).trim();

            ArgumentDefinition argumentDefinition = definitions.findArgumentDefinition(argument,definitionMatcher);
            embeddedValueMatch = new ArgumentMatch( argument, argumentDefinition, position );

            if( embeddedValueMatch != null && value != null )
                embeddedValueMatch.addValue( position, value );
        }

        // Prefer the standalone match...
        ArgumentMatch bestMatch = standaloneMatch;

        // ...But if the embedded value match is clearly better, choose it as the best match instead.
        if( (standaloneMatch == null || standaloneMatch.definition == null) &&
            (embeddedValueMatch != null && embeddedValueMatch.definition != null) )
            bestMatch = embeddedValueMatch;

        return bestMatch;
    }
}
