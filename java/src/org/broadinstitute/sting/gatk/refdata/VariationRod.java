package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.genotype.Variation;

/**
 * @author aaron
 *         <p/>
 *         Interface VariationRod
 *         <p/>
 *         This interface combines two interfaces: Variation and ReferenceOrderedDatum.  This
 * was required so that the reference ordered data require attribute would have an interface
 * that both specified variation and ROD compliance.   
 */
public interface VariationRod extends Variation, ReferenceOrderedDatum {
}
