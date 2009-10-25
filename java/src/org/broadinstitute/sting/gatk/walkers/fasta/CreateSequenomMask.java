package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.Iterator;

public class CreateSequenomMask extends RodWalker<Integer, Integer> {

    public void initialize() {}

	public Integer map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        int result = 0;
        if ( rodData == null ) // apparently, RodWalkers make funky map calls
            return 0;

        Iterator<ReferenceOrderedDatum> rods = rodData.getAllRods().iterator();
        while (rods.hasNext()) {
            ReferenceOrderedDatum rod = rods.next();
            if ( rod instanceof Variation && ((Variation)rod).isSNP() ) {
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