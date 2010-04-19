/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.commandline;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * Holds a pattern, along with how to get to the argument definitions that could match that pattern.
 */
public abstract class ParsingMethod {
    private final Pattern pattern;
    private final DefinitionMatcher definitionMatcher;

    /**
     * Create a new parsing method with the given identifying / validating pattern and definition matcher.
     * @param pattern The pattern
     * @param definitionMatcher The definition matcher.
     */
    private ParsingMethod( Pattern pattern, DefinitionMatcher definitionMatcher ) {
        this.pattern = pattern;
        this.definitionMatcher = definitionMatcher;
    }

    /**
     * Can the given token be parsed by this parsing method?
     * @param token Token to validate.
     * @return True if the given token matches.
     */
    public boolean matches( String token ) {
        Matcher matcher = pattern.matcher(token);
        return matcher.matches();        
    }

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

    public static ParsingMethod FullNameParsingMethod = new ParsingMethod(Pattern.compile("\\s*--([\\w\\.\\-]+)\\s*"),
                                                                          ArgumentDefinitions.FullNameDefinitionMatcher) {};
    public static ParsingMethod ShortNameParsingMethod = new ParsingMethod(Pattern.compile("\\s*-([\\w\\-]+)\\s*"),
                                                                           ArgumentDefinitions.ShortNameDefinitionMatcher) {};
}
