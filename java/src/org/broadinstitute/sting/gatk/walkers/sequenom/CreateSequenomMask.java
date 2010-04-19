package org.broadinstitute.sting.gatk.walkers.sequenom;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;

/**
 * Create a mask for use with the PickSequenomProbes walker.
 */
public class CreateSequenomMask extends RodWalker<Integer, Integer> {

    public void initialize() {}

	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        int result = 0;
        for ( VariantContext vc : tracker.getAllVariantContexts(ref) ) {
            if ( vc.isSNP() ) {
                out.println(context.getLocation());
                result = 1;
                break;
            }
        }

        return result;
    }

    public Integer reduceInit() {
        return 0;
    }

	public Integer reduce(Integer value, Integer sum) {
		return value + sum;
	}

    public void onTraversalDone(Integer sum) {
        logger.info("Found " + sum + " masking sites.");
    }
}