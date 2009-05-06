package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;

import java.lang.reflect.Field;
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

    static AliasProvider ShortNameAliasProvider = new AliasProvider() {
        /**
         * Short names can come in the form -Ofoo.txt, -O foo.txt, or -out (multi-character short name).
         * Given the argument name and built-in provided, see if these can be formed into some other argument
         * name.
         * @param argument Name of the argument, as parsed.  For a short name, will be a single letter.
         * @param value Value of the argument, as parsed.
         * @return Any potential aliases for the given shortname.
         */
        public List<String> getAliases( String argument, String value ) {
            List<String> aliases = new ArrayList<String>();
            aliases.add(argument+value);
            aliases.add(argument);
            return aliases;
        }

        /**
         * Is the value part of the given alias, or something separate that should be treated as an argument value.
         * @param alias The alias to use.
         * @param argument The parsed argument.
         * @param value The parsed value.
         * @return True if this alias should be used instead of the given value.
         */
        public boolean doesAliasConsumeValue( String alias, String argument, String value ) {
            return alias.equals(argument + value);
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
    public final Collection<ArgumentDefinition> argumentDefinitions;

    public ArgumentDefinitionGroup( String groupName, Collection<ArgumentDefinition> argumentDefinitions ) {
        this.groupName = groupName;
        this.argumentDefinitions = Collections.unmodifiableCollection( argumentDefinitions );
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
 * A specific argument definition.  Maps one-to-one with a field in some class.
 */
class ArgumentDefinition {
    /**
     * Full name of the argument.  Must have a value.
     */
    public final String fullName;

    /**
     * Short name of the argument.  Can be null.
     */
    public final String shortName;

    /**
     * Doc string for the argument.  Displayed in help.
     */
    public final String doc;

    /**
     * Is this argument required?
     */
    public final boolean required;

    /**
     * Is this argument exclusive of other arguments?
     */
    public final String exclusive;

    public final Class sourceClass;
    public final Field sourceField;

    /**
     * Creates a new argument definition.
     * @param argument Attributes of the argument, read from the source field.
     * @param sourceClass Source class for the argument, provided to the ParsingEngine.
     * @param sourceField Source field for the argument, extracted from the sourceClass.
     */
    public ArgumentDefinition( Argument argument, Class sourceClass, Field sourceField ) {
        this.sourceClass = sourceClass;
        this.sourceField = sourceField;

        fullName = argument.fullName().trim().length() > 0 ? argument.fullName().trim() : sourceField.getName().toLowerCase();
        shortName = argument.shortName().trim().length() > 0 ? argument.shortName().trim() : null;
        doc = argument.doc();
        required = argument.required() && !isFlag();
        exclusive = argument.exclusive().trim().length() > 0 ? argument.exclusive().trim() : null;
    }

    /**
     * Can this argument support multiple values, or just one?
     * @return True if the argument supports multiple values.
     */
    public boolean isMultiValued() {
        Class argumentType = sourceField.getType();
        return Collection.class.isAssignableFrom(argumentType) || sourceField.getType().isArray();         
    }

    /**
     * Is this argument a flag (meaning a boolean value whose presence indicates 'true').
     * @return True if this argument is a flag.
     */
    public boolean isFlag() {
        return (sourceField.getType() == Boolean.class) || (sourceField.getType() == Boolean.TYPE);
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

/**
 * A way to get alternate names for the argument given the recognized name and value.
 */
interface AliasProvider {
    /**
     * Give all alternate names for the given argument / value pair.  The aliases should
     * be returned in 'preferred order'.
     * @param argument The argument.
     * @param value The value.
     * @return All possible names.
     */
    List<String> getAliases( String argument, String value );

    /**
     * True if this alias 'consumes' the value, meaning that the argument + value together
     * represent some other alias.
     * @return True if the value should still be used.  False otherwise.
     */
    boolean doesAliasConsumeValue( String alias, String argument, String value );
}
