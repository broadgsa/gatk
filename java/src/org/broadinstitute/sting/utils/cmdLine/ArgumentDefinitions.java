package org.broadinstitute.sting.utils.cmdLine;

import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;

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
class ArgumentDefinitions {
    /**
     * Backing data set of argument stored by short name and long name.
     */
    private Map<String,ArgumentDefinition> argumentsByShortName = new HashMap<String,ArgumentDefinition>();
    private Map<String,ArgumentDefinition> argumentsByLongName = new HashMap<String,ArgumentDefinition>();

    /**
     * Returns the argument with the given short name.
     * @param shortName Argument short name.
     * @return The argument definition, or null if nothing matches.
     */
    public ArgumentDefinition getArgumentWithShortName( String shortName ) {
        return argumentsByShortName.get( shortName );
    }

    /**
     * Returns the argument with the given short name.
     * @param longName Argument long name.
     * @return The argument definition, or null if nothing matches.
     */
    public ArgumentDefinition getArgumentWithLongName( String longName ) {
        return argumentsByLongName.get( longName );
    }

    /**
     * Adds an argument to the this argument definition list.
     * @param argument The argument to add.
     * @param sourceClass Class where the argument was defined.
     * @param sourceField Field in which the argument was defined.
     */
    public void add( Argument argument, Class sourceClass, Field sourceField ) {
        ArgumentDefinition definition = new ArgumentDefinition( argument, sourceClass, sourceField );
        String fullName = argument.fullName().trim();
        String shortName = argument.shortName().trim();

        if( fullName.length() == 0 )
            throw new IllegalArgumentException( "Argument cannot have 0-length fullname." );

        argumentsByLongName.put( fullName, definition );
        if( shortName.length() != 0 )
            argumentsByShortName.put( shortName, definition );
    }

}

/**
 * A specific argument definition.  Maps one-to-one with a field in some class.
 */
class ArgumentDefinition {
    public final Argument argument;
    public final Class sourceClass;
    public final Field sourceField;

    /**
     * Creates a new argument definition.
     * @param argument Attributes of the argument, read from the source field.
     * @param sourceClass Source class for the argument, provided to the ParsingEngine.
     * @param sourceField Source field for the argument, extracted from the sourceClass.
     */
    public ArgumentDefinition( Argument argument, Class sourceClass, Field sourceField ) {
        this.argument = argument;
        this.sourceClass = sourceClass;
        this.sourceField = sourceField;
    }
}

/**
 * A general purpose accessor interface for ArgumentDefinitions.
 */
interface DefinitionMatcher {
    ArgumentDefinition get( ArgumentDefinitions argumentDefinitions, String key );
}

class FullNameDefinitionMatcher implements DefinitionMatcher {
    public ArgumentDefinition get( ArgumentDefinitions argumentDefinitions, String key ) {
        return argumentDefinitions.getArgumentWithLongName( key );
    }
}

class ShortNameDefinitionMatcher implements DefinitionMatcher {
    public ArgumentDefinition get( ArgumentDefinitions argumentDefinitions, String key ) {
        return argumentDefinitions.getArgumentWithShortName( key );
    }
}