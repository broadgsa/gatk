package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.utils.StingException;

import java.util.List;
import java.util.Collections;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * A group of argument definitions.
 */
public class ArgumentDefinitionGroup implements Iterable<ArgumentDefinition> {
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
