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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

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
        if( this.groupName == null )
            return other.groupName == null;
        return this.groupName.equals(other.groupName);
    }

    /**
     * Merges another argument group into this argument group.  Return a new
     * group since argument groups are supposed to be immutable. Asserts that
     * both argument groups have the same name.
     */
    public ArgumentDefinitionGroup merge( ArgumentDefinitionGroup other ) {
        if( !groupNameMatches(other) )
            throw new ReviewedGATKException("Unable to merge two argument groups with differing names.");

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

    /**
     * Reports whether all the arguments in this group are hidden.
     * @return True if all are hidden, false if some or none are hidden.
     */
    public boolean allHidden() {
        for(ArgumentDefinition argumentDefinition: argumentDefinitions) {
            if(!argumentDefinition.isHidden)
                return false;
        }
        return true;
    }
}
