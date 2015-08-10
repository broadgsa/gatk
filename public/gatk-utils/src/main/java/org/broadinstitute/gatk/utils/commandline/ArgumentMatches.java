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

import java.util.*;
/**
 * Represents a list of potential matches between the arguments defined
 * by the argument sources and the arguments passed in via the command line.
 */
public class ArgumentMatches implements Iterable<ArgumentMatch> {
    /**
     * Collection matches from argument definition to argument value.
     * Package protected access is deliberate.
     */
    Map<ArgumentMatchSite,ArgumentMatch> argumentMatches = new TreeMap<ArgumentMatchSite,ArgumentMatch>();

    /**
     * Provide a place to put command-line argument values that don't seem to belong to
     * any particular command-line option.
     */
    ArgumentMatch MissingArgument = new ArgumentMatch();

    /**
     * Get an iterator cycling through *unique* command-line argument <-> definition matches.
     * @return Iterator over all argument matches.
     */
    public Iterator<ArgumentMatch> iterator() {
        return getUniqueMatches().iterator();
    }

    /**
     * Create an empty ArgumentMatches object.
     */
    public ArgumentMatches() {
    }

    /**
     * Create a singleton ArgumentMatches object.
     * @param match Match to incorporate.
     */
    public ArgumentMatches( ArgumentMatch match ) {
        mergeInto( match );
    }

    /**
     * Returns the number of matches in this structure.
     * @return Count of the matches in this structure.
     */
    public int size() {
        return argumentMatches.size();
    }

    /**
     * Indicates whether the site contains a matched argument.
     * @param site Site at which to check.
     * @return True if the site has a match.  False otherwise.
     */
    boolean hasMatch( ArgumentMatchSite site ) {
        return argumentMatches.containsKey( site );
    }

    /**
     * Gets the match at a given site.
     * @param site Site at which to look for a match.
     * @return The match present at the given site.
     * @throws IllegalArgumentException if site does not contain a match.
     */
    ArgumentMatch getMatch( ArgumentMatchSite site ) {
        if( !argumentMatches.containsKey(site) )
            throw new IllegalArgumentException( "Site does not contain an argument: " + site );
        return argumentMatches.get(site);
    }

    /**
     * Does the match collection have a match for this argument definition.
     * @param definition Definition to match.
     * @return True if a match exists; false otherwise.
     */
    boolean hasMatch( ArgumentDefinition definition ) {
        return findMatches( definition ).size() > 0;
    }

    /**
     * Return all argument matches of this source.
     * @param parsingEngine Parsing engine.
     * @param argumentSource Argument source to match.
     * @return List of all matches.
     */

    ArgumentMatches findMatches(ParsingEngine parsingEngine, ArgumentSource argumentSource) {
        List<ArgumentDefinition> sourceDefinitions = parsingEngine.selectBestTypeDescriptor(argumentSource.field.getType()).createArgumentDefinitions(argumentSource);

        ArgumentMatches matches = new ArgumentMatches();
        for( ArgumentMatch argumentMatch: getUniqueMatches() ) {
            if( sourceDefinitions.contains(argumentMatch.definition) )
                matches.mergeInto( argumentMatch );
        }
        return matches;
    }

    /**
     * Return all argument matches of this definition.
     * @param definition Argument definition to match.
     * @return List of all matches.
     */
    ArgumentMatches findMatches( ArgumentDefinition definition ) {
        ArgumentMatches matches = new ArgumentMatches();
        for( ArgumentMatch argumentMatch: argumentMatches.values() ) {
            if( argumentMatch.definition == definition )
                matches.mergeInto( argumentMatch );
        }
        return matches;
    }

    /**
     * Find all successful matches (a 'successful' match is one paired with a definition).
     * @return All successful matches.
     */
    ArgumentMatches findSuccessfulMatches() {
        ArgumentMatches matches = new ArgumentMatches();
        for( ArgumentMatch argumentMatch: argumentMatches.values() ) {
            if( argumentMatch.definition != null )
                matches.mergeInto( argumentMatch );
        }
        return matches;
    }

    /**
     * Find arguments that are unmatched to any definition.
     * @return Set of matches that have no associated definition.
     */
    ArgumentMatches findUnmatched() {
        ArgumentMatches matches = new ArgumentMatches();
        for( ArgumentMatch argumentMatch: argumentMatches.values() ) {
            if( argumentMatch.definition == null )
                matches.mergeInto( argumentMatch );
        }
        return matches;
    }

    /**
     * Reformat the given entries with the given multiplexer and key.
     * TODO: Generify this.
     * @param multiplexer Multiplexer that controls the transformation process.
     * @param key Key which specifies the transform.
     * @return new argument matches.
     */
    ArgumentMatches transform(Multiplexer multiplexer, Object key) {
        ArgumentMatches newArgumentMatches = new ArgumentMatches();
        for(ArgumentMatch match: argumentMatches.values())
            newArgumentMatches.mergeInto(match.transform(multiplexer,key));
        return newArgumentMatches;
    }

    /**
     * Merges the given argument match into the set of existing argument matches.
     * If multiple arguments are present, those arguments will end up grouped.
     * @param match The match to merge into.
     */
    void mergeInto( ArgumentMatch match ) {
        boolean definitionExists = false;

        // Clone the list of argument matches to avoid ConcurrentModificationExceptions.
        for( ArgumentMatch argumentMatch: getUniqueMatches() ) {
            if( argumentMatch.definition == match.definition && argumentMatch.tags.equals(match.tags) ) {
                argumentMatch.mergeInto( match );
                for( ArgumentMatchSite site: match.sites.keySet() )
                    argumentMatches.put( site, argumentMatch );
                definitionExists = true;
            }
        }

        if( !definitionExists ) {
            for( ArgumentMatchSite site: match.sites.keySet() )
                argumentMatches.put( site, match );
        }
    }    

    /**
     * Determines, of the argument matches by position, which are unique and returns that list.
     * @return A unique set of matches.
     */
    private Set<ArgumentMatch> getUniqueMatches() {
        return new LinkedHashSet<ArgumentMatch>( argumentMatches.values() );
    }    
}
