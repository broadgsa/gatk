package org.broadinstitute.sting.utils.distributedutils;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 1/19/11
 * Time: 8:06 AM
 *
 * Information about processing locations and their owners.  Contains two basic data, associated
 * together.  The first is a genome loc, and the second is the name of the owner, as a string.
 *
 * chr1:1-10 Mark
 * chr2:11-20 DePristo
 *
 * would be two ProcessingLocs that first indicate that the first 10 bp of chr1 are owned by Mark,
 * and the second is owned by DePristo.
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
        this.owner = owner.intern();    // reduce memory consumption by interning the string
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    public String getOwner() {
        return owner;
    }

    /**
     * Returns true iff the owner of this processing loc is name.  Can be used to determine
     * the owner of this processing location.
     *
     * @param name
     * @return
     */
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
