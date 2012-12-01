package org.broadinstitute.sting.utils.activeregion;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 11/26/12
 * Time: 2:35 PM
 *
 * Describes how a read relates to an assigned ActiveRegion
 */
public enum ActiveRegionReadState {
    PRIMARY,        // This is the read's primary region
    NONPRIMARY,     // This region overlaps the read, but it is not primary
    EXTENDED,       // This region would overlap the read if it were extended
    UNMAPPED        // This read is not mapped
}
