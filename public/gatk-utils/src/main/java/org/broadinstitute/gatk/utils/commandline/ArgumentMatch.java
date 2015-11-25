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
 * A mapping of all the sites where an argument definition maps to a site on the command line.
 */
public class ArgumentMatch implements Iterable<ArgumentMatch> {
    /**
     * The argument definition that's been matched.
     */
    public final ArgumentDefinition definition;

    /**
     * The text that's been matched, as it appears in the command line arguments.
     */
    public final String label;

    /**
     * Maps indices of command line arguments to values paired with that argument.
     */
    public final SortedMap<ArgumentMatchSite,List<ArgumentMatchValue>> sites = new TreeMap<ArgumentMatchSite,List<ArgumentMatchValue>>();

    /**
     * An ordered, freeform collection of tags.
     */
    public final Tags tags;

    /**
     * Create a new argument match, defining its properties later.  Used to create invalid arguments.
     */
    public ArgumentMatch() {
        this(null,null);
    }

    /**
     * Minimal constructor for transform function.
     * @param label Label of the argument match.  Must not be null.
     * @param definition The associated definition, if one exists.  May be null.
     */
    private ArgumentMatch(final String label, final ArgumentDefinition definition) {
        this.label = label;
        this.definition = definition;
        this.tags = new Tags();
    }

    /**
     * A simple way of indicating that an argument with the given label and definition exists at this site.
     * @param label Label of the argument match.  Must not be null.
     * @param definition The associated definition, if one exists.  May be null.
     * @param site Position of the argument.  Must not be null.
     * @param tags ordered freeform text tags associated with this argument.
     */
    public ArgumentMatch(final String label, final ArgumentDefinition definition, final ArgumentMatchSite site, final Tags tags) {
        this( label, definition, site, null, tags );
    }

    /**
     * A simple way of indicating that an argument with the given label and definition exists at this site.
     * @param label Label of the argument match.  Must not be null.
     * @param definition The associated definition, if one exists.  May be null.
     * @param site Position of the argument.  Must not be null.
     * @param value Value for the argument at this position.
     * @param tags ordered freeform text tags associated with this argument.
     */
    private ArgumentMatch(final String label, final ArgumentDefinition definition, final ArgumentMatchSite site, final ArgumentMatchValue value, final Tags tags) {
        this.label = label;
        this.definition = definition;

        ArrayList<ArgumentMatchValue> values = new ArrayList<ArgumentMatchValue>();
        if( value != null )
            values.add(value);
        sites.put(site,values );

        this.tags = tags;
    }

    /**
     * Check to see whether two ArgumentMatch objects are equal.
     * @param other Object to which this should be compared.
     * @return True if objects are equal, false if objects are not equal or incomparable.
     */
    @Override
    public boolean equals(Object other) {
        // this clearly isn't null, since this.equals() when this == null would result in an NPE.
        if(other == null)
            return false;
        if(!(other instanceof ArgumentMatch))
            return false;
        ArgumentMatch otherArgumentMatch = (ArgumentMatch)other;
        return this.definition.equals(otherArgumentMatch.definition) &&
                this.label.equals(otherArgumentMatch.label) &&
                this.sites.equals(otherArgumentMatch.sites) &&
                this.tags.equals(otherArgumentMatch.tags);
    }


    /**
     * Reformat the given entries with the given multiplexer and key.
     * TODO: Generify this.
     * @param multiplexer Multiplexer that controls the transformation process.
     * @param key Key which specifies the transform.
     * @return A variant of this ArgumentMatch with all keys transformed.
     */
    @SuppressWarnings("unchecked")
    ArgumentMatch transform(Multiplexer multiplexer, Object key) {
        SortedMap<ArgumentMatchSite,List<ArgumentMatchValue>> newIndices = new TreeMap<ArgumentMatchSite,List<ArgumentMatchValue>>();
        for(Map.Entry<ArgumentMatchSite,List<ArgumentMatchValue>> site: sites.entrySet()) {
            List<ArgumentMatchValue> newEntries = new ArrayList<ArgumentMatchValue>();
            for(ArgumentMatchValue entry: site.getValue())
                newEntries.add(new ArgumentMatchStringValue(multiplexer.transformArgument(key,entry.asString())));
            newIndices.put(site.getKey(),newEntries);
        }
        ArgumentMatch newArgumentMatch = new ArgumentMatch(label,definition);
        newArgumentMatch.sites.putAll(newIndices);
        return newArgumentMatch;
    }

    /**
     * Return a string representation of the given argument match, for debugging purposes.
     * @return String representation of the match.
     */
    public String toString() {
        return label;
    }

    /**
     * Creates an iterator that walks over each individual match at each position of a given argument.
     * @return An iterator over the individual matches in this argument.  Will not be null.
     */
    public Iterator<ArgumentMatch> iterator() {
        return new Iterator<ArgumentMatch>() {
            /**
             * Iterate over each the available site.
             */
            private Iterator<ArgumentMatchSite> siteIterator = null;

            /**
             * Iterate over each available token.
             */
            private Iterator<ArgumentMatchValue> tokenIterator = null;

            /**
             * The next site to return.  Null if none remain.
             */
            ArgumentMatchSite nextSite = null;

            /**
             * The next token to return.  Null if none remain.
             */
            ArgumentMatchValue nextToken = null;

            {
                siteIterator = sites.keySet().iterator();
                prepareNext();
            }

            /**
             * Is there a nextToken available to return?
             * @return True if there's another token waiting in the wings.  False otherwise.
             */
            public boolean hasNext() {
                return nextSite != null;
            }

            /**
             * Get the next token, if one exists.  If not, throw an IllegalStateException.
             * @return The next ArgumentMatch in the series.  Should never be null.
             */
            public ArgumentMatch next() {
                if( nextSite == null )
                    throw new IllegalStateException( "No more ArgumentMatches are available" );

                ArgumentMatch match = new ArgumentMatch( label, definition, nextSite, nextToken, tags );
                prepareNext();
                return match;
            }

            /**
             * Initialize the next ArgumentMatch to return.  If no ArgumentMatches are available,
             * initialize nextSite / nextToken to null.
             */
            private void prepareNext() {
                if( tokenIterator != null && tokenIterator.hasNext() ) {
                    nextToken = tokenIterator.next();
                }
                else {
                    nextSite = null;
                    nextToken = null;

                    // Do a nested loop.  While more data is present in the inner loop, grab that data.
                    // Otherwise, troll the outer iterator looking for more data.
                    while( siteIterator.hasNext() ) {
                        nextSite = siteIterator.next();
                        if( sites.get(nextSite) != null ) {
                            tokenIterator = sites.get(nextSite).iterator();
                            nextToken = tokenIterator.hasNext() ? tokenIterator.next() : null;
                            break;
                        }
                    }
                }

            }

            /**
             * Remove is unsupported in this context.
             */
            public void remove() {
                throw new UnsupportedOperationException("Cannot remove an argument match from the collection while iterating.");
            }
        };
    }

    /**
     * Merge two ArgumentMatches, so that the values for all arguments go into the
     * same data structure.
     * @param other The other match to merge into.
     */
    public void mergeInto( ArgumentMatch other ) {
        sites.putAll(other.sites);
    }

    /**
     * Associate a value with this merge maapping.
     * @param site site of the command-line argument to which this value is mated.
     * @param value Text representation of value to add.
     */
    public void addValue( ArgumentMatchSite site, ArgumentMatchValue value ) {
        if( !sites.containsKey(site) || sites.get(site) == null )
            sites.put(site, new ArrayList<ArgumentMatchValue>() );
        sites.get(site).add(value);
    }

    /**
     * Does this argument already have a value at the given site?
     * Arguments are only allowed to be single-valued per site, and
     * flags aren't allowed a value at all.
     * @param site Site at which to check for values.
     * @return True if the argument has a value at the given site.  False otherwise.
     */
    public boolean hasValueAtSite( ArgumentMatchSite site ) {
        return (sites.get(site) != null && sites.get(site).size() >= 1) || isArgumentFlag();
    }

    /**
     * Return the values associated with this argument match.
     * @return A collection of the string representation of these value.
     */
    public List<ArgumentMatchValue> values() {
        final List<ArgumentMatchValue> values = new ArrayList<ArgumentMatchValue>();
        for ( final List<ArgumentMatchValue> siteValue : sites.values() ) {
            if ( siteValue != null )
                values.addAll(siteValue);
            else
                values.add(null);
        }
        return values;
    }

    /**
     * Convenience method returning true if the definition is a flag.
     * @return True if definition is known to be a flag; false if not known to be a flag.
     */
    private boolean isArgumentFlag() {
        return definition != null && definition.isFlag;
    }
}
