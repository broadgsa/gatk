package org.broadinstitute.sting.utils.threading;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 1/19/11
 * Time: 8:06 AM
 *
 * Information about processing locations and their owners
 */
public class ProcessingLoc implements HasGenomeLocation {
    private final GenomeLoc loc;
    private final String owner;

    /**
     * Create a loc that's already owned
     * @param loc
     * @param owner
     */
    public ProcessingLoc(GenomeLoc loc, String owner) {
        if ( loc == null || owner == null ) {
            throw new ReviewedStingException("BUG: invalid ProcessingLoc detected: " + loc + " owner " + owner);
        }

        this.loc = loc;
        this.owner = owner;
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    public String getOwner() {
        return owner;
    }

    public boolean isOwnedBy(String name) {
        return getOwner().equals(name);
    }

    public String toString() { return String.format("ProcessingLoc(%s,%s)", loc, owner); }

    public boolean equals(Object other) {
        if (other instanceof ProcessingLoc )
            return this.loc.equals(((ProcessingLoc)other).loc) && this.owner.equals(((ProcessingLoc)other).owner);
        else
            return false;
    }

    public int compareTo(ProcessingLoc other) {
        return this.getLocation().compareTo(other.getLocation());
    }
}
