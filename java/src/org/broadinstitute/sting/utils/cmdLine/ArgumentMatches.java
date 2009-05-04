package org.broadinstitute.sting.utils.cmdLine;

import java.util.Iterator;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: May 3, 2009
 * Time: 6:36:43 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */

/**
 * Represents a list of potential matches between the arguments defined
 * by the argument sources and the arguments passed in via the command line.
 */
public class ArgumentMatches implements Iterable<ArgumentMatch> {
    /**
     * Collection matches from argument definition to argument value.
     * Package protected access is deliberate.
     */
    Map<Integer,ArgumentMatch> argumentMatches = new TreeMap<Integer,ArgumentMatch>();

    void mergeInto( ArgumentMatch match ) {
        boolean definitionExists = false;

        // Clone the list of argument matches to avoid ConcurrentModificationExceptions.
        Set<ArgumentMatch> uniqueMatches = getUniqueMatches();
        for( ArgumentMatch argumentMatch: uniqueMatches ) {
            if( argumentMatch.definition.equals(match.definition) ) {
                argumentMatch.mergeInto( match );
                for( int index: match.indices )
                    argumentMatches.put( index, argumentMatch );
                definitionExists = true;
            }
        }

        if( !definitionExists ) {
            for( int index: match.indices )
                argumentMatches.put( index, match );
        }
    }

    /**
     * Get an iterator cycling through *unique* command-line argument <-> definition matches.
     * @return Iterator over all argument matches.
     */
    public Iterator<ArgumentMatch> iterator() {
        return getUniqueMatches().iterator();
    }

    /**
     * Indicates whether the site contains a matched argument.
     * @param site Site at which to check.
     * @return True if the site has a match.  False otherwise.
     */
    public boolean hasMatch( int site ) {
        return argumentMatches.containsKey( site );
    }

    /**
     * Gets the match at a given site.
     * @param site Site at which to look for a match.
     * @return The match present at the given site.
     * @throws IllegalArgumentException if site does not contain a match.
     */
    public ArgumentMatch getMatch( int site ) {
        if( !argumentMatches.containsKey(site) )
            throw new IllegalArgumentException( "Site does not contain an argument: " + site );
        return argumentMatches.get(site);
    }

    /**
     * Determines, of the argument matches by position, which are unique and returns that list.
     * @return A unique set of matches.
     */
    private Set<ArgumentMatch> getUniqueMatches() {
        return new HashSet<ArgumentMatch>( argumentMatches.values() );
    }
}

/**
 * A mapping of all the sites where an argument definition maps to a site on the command line.
 */
class ArgumentMatch {
    /**
     * The argument definition that's been matched.
     */
    public final ArgumentDefinition definition;

    /**
     * Index into the string of arguments where this match was found.
     */
    public final Set<Integer> indices = new HashSet<Integer>();

    /**
     * The values associated with this parameter.
     */
    public final List<String> values = new ArrayList<String>();

    public ArgumentMatch( ArgumentDefinition definition, int index ) {
        this.definition = definition;
        indices.add(index);
    }

    /**
     * Merge two ArgumentMatches, so that the values for all arguments go into the
     * same data structure.
     * @param other The other match to merge into.
     */
    public void mergeInto( ArgumentMatch other ) {
        indices.addAll(other.indices);
        values.addAll(other.values);
    }

    /**
     * Associate a value with this merge maapping.
     * @param value Text representation of value to add.
     */
    public void addValue( String value ) {
        this.values.add(value);
    }
}