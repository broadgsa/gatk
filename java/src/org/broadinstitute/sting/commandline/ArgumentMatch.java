/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.commandline;

import org.broadinstitute.sting.gatk.walkers.Multiplexer;

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
     * Maps indicies of command line arguments to values paired with that argument.
     */
    public final SortedMap<Integer,List<String>> indices = new TreeMap<Integer,List<String>>();

    /**
     * An ordered, freeform collection of tags.
     */
    public final List<String> tags;

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
    private ArgumentMatch(String label,ArgumentDefinition definition) {
        this.label = label;
        this.definition = definition;
        this.tags = Collections.emptyList();
    }

    /**
     * A simple way of indicating that an argument with the given label and definition exists at this index.
     * @param label Label of the argument match.  Must not be null.
     * @param definition The associated definition, if one exists.  May be null.
     * @param index Position of the argument.  Must not be null.
     * @param tags ordered freeform text tags associated with this argument.
     */
    public ArgumentMatch( String label, ArgumentDefinition definition, int index, List<String> tags ) {
        this( label, definition, index, null, tags );
    }

    /**
     * A simple way of indicating that an argument with the given label and definition exists at this index.
     * @param label Label of the argument match.  Must not be null.
     * @param definition The associated definition, if one exists.  May be null.
     * @param index Position of the argument.  Must not be null.
     * @param value Value for the argument at this position.
     * @param tags ordered freeform text tags associated with this argument.
     */
    private ArgumentMatch( String label, ArgumentDefinition definition, int index, String value, List<String> tags ) {
        this.label = label;
        this.definition = definition;

        ArrayList<String> values = new ArrayList<String>();
        if( value != null )
            values.add(value);
        indices.put(index,values );

        this.tags = tags;
    }

    /**
     * Reformat the given entries with the given multiplexer and key.
     * TODO: Generify this.
     * @param multiplexer Multiplexer that controls the transformation process.
     * @param key Key which specifies the transform.
     * @return A variant of this ArgumentMatch with all keys transformed.
     */
    ArgumentMatch transform(Multiplexer multiplexer, Object key) {
        SortedMap<Integer,List<String>> newIndices = new TreeMap<Integer,List<String>>();
        for(Map.Entry<Integer,List<String>> index: indices.entrySet()) {
            List<String> newEntries = new ArrayList<String>();
            for(String entry: index.getValue())
                newEntries.add(multiplexer.transformArgument(key,entry));
            newIndices.put(index.getKey(),newEntries);
        }
        ArgumentMatch newArgumentMatch = new ArgumentMatch(label,definition);
        newArgumentMatch.indices.putAll(newIndices);
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
             * Iterate over each the available index.
             */
            private Iterator<Integer> indexIterator = null;

            /**
             * Iterate over each available token.
             */
            private Iterator<String> tokenIterator = null;

            /**
             * The next index to return.  Null if none remain.
             */
            Integer nextIndex = null;

            /**
             * The next token to return.  Null if none remain.
             */
            String nextToken = null;

            {
                indexIterator = indices.keySet().iterator();
                prepareNext();
            }

            /**
             * Is there a nextToken available to return?
             * @return True if there's another token waiting in the wings.  False otherwise.
             */
            public boolean hasNext() {
                return nextToken != null;    
            }

            /**
             * Get the next token, if one exists.  If not, throw an IllegalStateException.
             * @return The next ArgumentMatch in the series.  Should never be null.
             */
            public ArgumentMatch next() {
                if( nextIndex == null || nextToken == null )
                    throw new IllegalStateException( "No more ArgumentMatches are available" );

                ArgumentMatch match = new ArgumentMatch( label, definition, nextIndex, nextToken, tags );
                prepareNext();
                return match;
            }

            /**
             * Initialize the next ArgumentMatch to return.  If no ArgumentMatches are available,
             * initialize nextIndex / nextToken to null.
             */
            private void prepareNext() {
                if( tokenIterator != null && tokenIterator.hasNext() ) {
                    nextToken = tokenIterator.next();
                }
                else {
                    nextIndex = null;
                    nextToken = null;

                    // Do a nested loop.  While more data is present in the inner loop, grab that data.
                    // Otherwise, troll the outer iterator looking for more data.
                    while( indexIterator.hasNext() ) {
                        nextIndex = indexIterator.next();
                        if( indices.get(nextIndex) != null ) {
                            tokenIterator = indices.get(nextIndex).iterator();
                            if( tokenIterator.hasNext() ) {
                                nextToken = tokenIterator.next();
                                break;
                            }
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
        return definition != null && definition.isFlag;
    }
}
