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

import org.broadinstitute.gatk.utils.Utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
    public ArgumentMatch match( ArgumentDefinitions definitions, String token, ArgumentMatchSite position ) {
        // If the argument is valid, parse out the argument.
        Matcher matcher = pattern.matcher(token);

        // Didn't match?  Must be bad input.
        if( !matcher.matches() )
            throw new IllegalArgumentException( String.format("Unable to parse token %s with pattern %s", token, pattern.pattern()) );

        String argument = matcher.group(1).trim();

        Tags tags = parseTags(argument, matcher.group(2));

        // Find the most appropriate argument definition for the given argument.
        ArgumentDefinition argumentDefinition = definitions.findArgumentDefinition( argument, definitionMatcher );

        // Try to find a matching argument.  If found, label that as the match.  If not found, add the argument
        // with a null definition.
        return new ArgumentMatch(argument,argumentDefinition,position,tags);
    }

    public static Tags parseTags(String argument, String tagString) {
        Tags tags = new Tags();
        if (tagString != null) {
            for(String tag: Utils.split(tagString, ",")) {
                // Check for presence of an '=' sign, indicating a key-value pair in the tag line.
                int equalDelimiterPos = tag.indexOf('=');
                if(equalDelimiterPos >= 0) {
                    // Sanity check; ensure that there aren't multiple '=' in this key-value pair.
                    if(tag.indexOf('=',equalDelimiterPos+1) >= 0)
                        throw new ArgumentException(String.format("Tag %s passed to argument %s is malformed.  Please ensure that " +
                                "key-value tags are of the form <key>=<value>, and neither key " +
                                "nor value contain the '=' character", tag, argument));
                    tags.addKeyValueTag(tag.substring(0,equalDelimiterPos),tag.substring(equalDelimiterPos+1));
                }
                else
                    tags.addPositionalTag(tag);

            }
        }
        return tags;
    }

    /**
     * A command-line argument always starts with an alphabetical character or underscore followed by any word character.
     */
    private static final String ARGUMENT_TEXT = "[A-Za-z_][\\w\\-\\.]*";

    /**
     * Tags, on the other hand, can start with any word character.
     */
    private static final String TAG_TEXT = "[\\w\\-\\.\\=]*";

    public static final ParsingMethod FullNameParsingMethod = new ParsingMethod(Pattern.compile(String.format("\\s*--(%1$s)(?:\\:(%2$s(?:,%2$s)*))?\\s*",ARGUMENT_TEXT,TAG_TEXT)),
                                                                          ArgumentDefinitions.FullNameDefinitionMatcher) {};
    public static final ParsingMethod ShortNameParsingMethod = new ParsingMethod(Pattern.compile(String.format("\\s*-(%1$s)(?:\\:(%2$s(?:,%2$s)*))?\\s*",ARGUMENT_TEXT,TAG_TEXT)),
                                                                           ArgumentDefinitions.ShortNameDefinitionMatcher) {};
}
