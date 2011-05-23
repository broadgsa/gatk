package org.broadinstitute.sting.utils;

import com.google.java.contract.Ensures;

/**
 * Indicates that this object has a genomic location and provides a systematic interface to get it.
 */
public interface HasGenomeLocation {
    @Ensures("result != null")
    public GenomeLoc getLocation();
}
