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

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * A collection of argument definitions.
 */
public class ArgumentDefinitions implements Iterable<ArgumentDefinition> {
    /**
     * Backing data set of argument stored by short name and long name.
     */
    private Set<ArgumentDefinition> argumentDefinitions = new HashSet<ArgumentDefinition>();

    /**
     * The groupings of argument definitions.  Used mainly for help output.
     */
    private Set<ArgumentDefinitionGroup> argumentDefinitionGroups = new HashSet<ArgumentDefinitionGroup>();

    /**
     * Adds an argument to the this argument definition list.
     * @param argumentDefinitionGroup The group of arguments to add.
     */
    public void add( ArgumentDefinitionGroup argumentDefinitionGroup ) {
        for( ArgumentDefinition definition: argumentDefinitionGroup ) {
            // Do some basic validation before adding the definition. 
            if( definition.fullName.length() == 0 )
                throw new IllegalArgumentException( "Argument cannot have 0-length fullname." );
            if( hasArgumentDefinition( definition.fullName, FullNameDefinitionMatcher ) )
                throw new ReviewedGATKException("Duplicate definition of argument with full name: " + definition.fullName);
            if( definition.shortName != null && hasArgumentDefinition( definition.shortName, ShortNameDefinitionMatcher ) )
                throw new ReviewedGATKException("Duplicate definition of argument with short name: " + definition.shortName);

            argumentDefinitions.add( definition );
        }

        // Find an existing argument definition group with this name.
        // If one exists, merge this group into the other.
        Iterator<ArgumentDefinitionGroup> definitionGroupIterator = argumentDefinitionGroups.iterator();
        while( definitionGroupIterator.hasNext() ) {
            ArgumentDefinitionGroup candidate = definitionGroupIterator.next();            
            if( candidate.groupNameMatches(argumentDefinitionGroup) ) {
                argumentDefinitionGroup = candidate.merge(argumentDefinitionGroup);
                definitionGroupIterator.remove();
            }
        }

        // Otherwise, add the new group.
        argumentDefinitionGroups.add( argumentDefinitionGroup );
    }

    /**
     * Are there any argument definitions matching the given property?
     * @param property Property to find.
     * @param matcher Method of matching a given property.
     * @return True if one or multiple argument definitions match; false otherwise.
     */
    public boolean hasArgumentDefinition( Object property, DefinitionMatcher matcher ) {
        return findArgumentDefinitions( property, matcher ).size() > 0;
    }

    /**
     * Find the given definition matching this property.
     * @param property Property to find.
     * @param matcher Method of matching a given property.
     * @return The ArgumentDefinition matching the given property.  Null if none matches.
     * @throws IllegalArgumentException if multiple arguments match this definition.
     */
    public ArgumentDefinition findArgumentDefinition( Object property, DefinitionMatcher matcher ) {
        Collection<ArgumentDefinition> selectedDefinitions = findArgumentDefinitions( property, matcher );
        if( selectedDefinitions.size() > 1 )
            throw new IllegalArgumentException("Multiple argument definitions match the selected property: " + property);

        if( selectedDefinitions.size() == 0 )
            return null;

        return selectedDefinitions.iterator().next();
    }

    /**
     * Find all argument definitions matching a certain category.
     * @param property Property to inspect.
     * @param matcher Test to see whether property matches.
     * @return All argument definitions matching a certain object.
     */
    public Collection<ArgumentDefinition> findArgumentDefinitions( Object property, DefinitionMatcher matcher ) {
        Set<ArgumentDefinition> selectedArgumentDefinitions = new HashSet<ArgumentDefinition>();
        for( ArgumentDefinition argumentDefinition: argumentDefinitions ) {
            if( matcher.matches( argumentDefinition, property ) )
                selectedArgumentDefinitions.add( argumentDefinition );
        }
        return selectedArgumentDefinitions;
    }

    /**
     * Return a list of the available argument groups.
     * @return All the argument groups that have been added.
     */
    public Collection<ArgumentDefinitionGroup> getArgumentDefinitionGroups() {
        return argumentDefinitionGroups;
    }

    /**
     * Iterates through all command-line arguments.
     * @return an iterator over command-line arguments.
     */
    public Iterator<ArgumentDefinition> iterator() {
        return argumentDefinitions.iterator();
    }

    /**
     * Match the full name of a definition.
     */
    static DefinitionMatcher FullNameDefinitionMatcher = new DefinitionMatcher() {
        public boolean matches( ArgumentDefinition definition, Object key ) {
            if( definition.fullName == null )
                return key == null;
            else
                return definition.fullName.equals( key );
        }        
    };

    /**
     * Match the short name of a definition.
     */
    static DefinitionMatcher ShortNameDefinitionMatcher = new DefinitionMatcher() {
        public boolean matches( ArgumentDefinition definition, Object key ) {
            if( definition.shortName == null )
                return key == null;
            else
                return definition.shortName.equals( key );
        }
    };

    /**
     * Find all required definitions.
     */
    static DefinitionMatcher RequiredDefinitionMatcher = new DefinitionMatcher() {
        public boolean matches( ArgumentDefinition definition, Object key ) {
            if( !(key instanceof Boolean) )
                throw new IllegalArgumentException("RequiredDefinitionMatcher requires boolean key");
            return definition.required == (Boolean)key;
        }
    };

    static DefinitionMatcher VerifiableDefinitionMatcher = new DefinitionMatcher() {
        public boolean matches( ArgumentDefinition definition, Object key ) {
            // We can perform some sort of validation for anything that isn't a flag or enum.
            // Because enums can have a default value, it might be valid to specify an enum argument with no value
            return !definition.isFlag  && !definition.argumentType.isEnum();
        }        
    };
}

/**
 * A Comparator-esque interface for finding argument definitions within a collection.
 */
interface DefinitionMatcher {
    /**
     * Does the given definition match the provided key?
     * @param definition The definition to inspect.
     * @param key The value to match.
     * @return True if the key matches the definition, false otherwise.
     */
    boolean matches( ArgumentDefinition definition, Object key );
}
