package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;

/**
 * PrintRODsWalker prints out all of the RODs that it sees (using the ROD's toString method)
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=ReferenceOrderedDatum.class))
public class PrintRODsWalker extends RefWalker<Integer, Integer> {

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    /**
     * For each site of interest, rescore the genotype likelihoods by applying the specified feature set.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ReferenceOrderedDatum variant = tracker.lookup("variant", null);

        if (variant != null )
            out.println(variant);

        return 1;
    }

    /**
     * Increment the number of rods processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of rods processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {}    
}