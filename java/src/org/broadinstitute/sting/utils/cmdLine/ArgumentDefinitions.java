package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;

import java.util.Set;
import java.util.HashSet;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Collections;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: May 3, 2009
 * Time: 6:02:04 PM
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
 * A collection of argument definitions.
 */
class ArgumentDefinitions implements Iterable<ArgumentDefinition> {
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
                throw new StingException("Duplicate definition of argument with full name: " + definition.fullName);
            if( definition.shortName != null && hasArgumentDefinition( definition.shortName, ShortNameDefinitionMatcher ) )
                throw new StingException("Duplicate definition of argument with short name: " + definition.shortName);

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
    Collection<ArgumentDefinitionGroup> getArgumentDefinitionGroups() {
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
            return definition.validation != null;
        }        
    };
}

/**
 * A group of argument definitions.
 */
class ArgumentDefinitionGroup implements Iterable<ArgumentDefinition> {
    /**
     * Name of this group.
     */
    public final String groupName;

    /**
     * The argument definitions associated with this group.
     */
    public final List<ArgumentDefinition> argumentDefinitions;

    public ArgumentDefinitionGroup( String groupName, List<ArgumentDefinition> argumentDefinitions ) {
        this.groupName = groupName;
        this.argumentDefinitions = Collections.unmodifiableList( argumentDefinitions );
    }

    /**
     * Does the name of this argument group match the name of another?
     */
    public boolean groupNameMatches( ArgumentDefinitionGroup other ) {
        if( this.groupName == null && other.groupName == null )
            return true;
        if( this.groupName == null && other.groupName != null )
            return false;
        return this.groupName.equals(other.groupName);
    }

    /**
     * Merges another argument group into this argument group.  Return a new
     * group since argument groups are supposed to be immutable. Asserts that
     * both argument groups have the same name.
     */
    public ArgumentDefinitionGroup merge( ArgumentDefinitionGroup other ) {
        if( !groupNameMatches(other) )
            throw new StingException("Unable to merge two argument groups with differing names.");

        // Create a merged definition group.
        List<ArgumentDefinition> mergedDefinitions = new ArrayList<ArgumentDefinition>();
        mergedDefinitions.addAll(this.argumentDefinitions);
        mergedDefinitions.addAll(other.argumentDefinitions);

        return new ArgumentDefinitionGroup(groupName,mergedDefinitions);
    }

    /**
     * Iterate over the arguments in an argument definition group.
     * @return
     */
    public Iterator<ArgumentDefinition> iterator() {
        return argumentDefinitions.iterator();
    }
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
