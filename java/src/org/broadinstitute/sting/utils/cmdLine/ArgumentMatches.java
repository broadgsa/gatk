package org.broadinstitute.sting.utils.cmdLine;

import java.util.Iterator;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.Collection;
import java.util.HashMap;
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

    /**
     * Provide a place to put command-line argument values that don't seem to belong to
     * any particular command-line option.
     */
    public ArgumentMatch MissingArgument = new ArgumentMatch();    

    void mergeInto( ArgumentMatch match ) {
        boolean definitionExists = false;

        // Clone the list of argument matches to avoid ConcurrentModificationExceptions.
        for( ArgumentMatch argumentMatch: getUniqueMatches() ) {
            if( argumentMatch.definition == match.definition ) {
                argumentMatch.mergeInto( match );
                for( int index: match.indices.keySet() )
                    argumentMatches.put( index, argumentMatch );
                definitionExists = true;
            }
        }

        if( !definitionExists ) {
            for( int index: match.indices.keySet() )
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

    /**
     * Does the match collection have a match for this argument definition.
     * @param definition Definition to match.
     * @return True if a match exists; false otherwise.
     */
    public boolean hasMatch( ArgumentDefinition definition ) {
        return findMatches( definition ).size() > 0;
    }

    /**
     * Return all argument matches of this definition.
     * @param definition Argument definition to match.
     * @return List of all matches.
     */
    public Collection<ArgumentMatch> findMatches( ArgumentDefinition definition ) {
        Collection<ArgumentMatch> matches = new HashSet<ArgumentMatch>();
        for( ArgumentMatch argumentMatch: getUniqueMatches() ) {
            if( argumentMatch.definition == definition )
                matches.add( argumentMatch );
        }        
        return matches;
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
     * The text that's been matched, as it appears in the command line arguments.
     */
    public final String label;

    /**
     * Maps indicies of command line arguments to values paired with that argument.
     */
    public final Map<Integer,List<String>> indices = new HashMap<Integer,List<String>>();

    /**
     * Create a new argument match, defining its properties later.  Used to create invalid arguments.
     */
    public ArgumentMatch() {
        definition = null;
        label = null;
    }

    public ArgumentMatch( String label, ArgumentDefinition definition, int index ) {
        this.label = label;
        this.definition = definition;
        indices.put(index,null);
    }

    /**
     * Merge two ArgumentMatches, so that the values for all arguments go into the
     * same data structure.
     * @param other The other match to merge into.
     */
    public void mergeInto( ArgumentMatch other ) {
        indices.putAll(other.indices);
    }

    /**
     * Associate a value with this merge maapping.
     * @param index index of the command-line argument to which this value is mated.
     * @param value Text representation of value to add.
     */
    public void addValue( int index, String value ) {
        if( !indices.containsKey(index) || indices.get(index) == null )
            indices.put(index, new ArrayList<String>() );
        indices.get(index).add(value);
    }

    /**
     * Does this argument already have a value at the given site?
     * Arguments are only allowed to be single-valued per site, and
     * flags aren't allowed a value at all.
     * @param index Index at which to check for values.
     * @return True if the argument has a value at the given site.  False otherwise.
     */
    public boolean hasValueAtSite( int index ) {
        return (indices.get(index) != null && indices.get(index).size() >= 1) || isArgumentFlag();
    }

    /**
     * Return the values associated with this argument match.
     * @return A collection of the string representation of these value.
     */
    public List<String> values() {
        List<String> values = new ArrayList<String>();
        for( int index: indices.keySet() ) {
            if( indices.get(index) != null )
                values.addAll(indices.get(index));
        }
        return values;
    }

    /**
     * Convenience method returning true if the definition is a flag.
     * @return True if definition is known to be a flag; false if not known to be a flag.
     */
    private boolean isArgumentFlag() {
        return definition != null && definition.isFlag();    
    }
}