package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

/**
 * @author aaron
 *         <p/>
 *         Interface Variant
 *         <p/>
 *         This class represents a variant
 */
public interface Variation {
    // the types of variants we currently allow
    public enum VARIANT_TYPE {
        SNP, INSERTION, DELETION, REFERENCE // though reference is not really a variant, we need to represent it
    }
}
