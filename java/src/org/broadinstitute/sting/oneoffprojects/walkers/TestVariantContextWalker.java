package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.EnumSet;

/**
 * Test routine for new VariantContext object
 */
public class TestVariantContextWalker extends RodWalker<Integer, Integer> {
    @Argument(fullName="takeFirstOnly", doc="Only take the first second at a locus, as opposed to all", required=false)
    boolean takeFirstOnly = false;

    @Argument(fullName="onlyContextsOfType", doc="Only take variant contexts of this type", required=false)
    VariantContext.Type onlyOfThisType = null;

    @Argument(fullName="onlyContextsStartinAtCurrentPosition", doc="Only take variant contexts at actually start at the current position, excluding those at span to the current location but start earlier", required=false)
    boolean onlyContextsStartinAtCurrentPosition = false;

    @Argument(fullName="printPerLocus", doc="If true, we'll print the variant contexts, in addition to counts", required=false)
    boolean printContexts = false;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( ref == null )
            return 0;
        else {
            EnumSet<VariantContext.Type> allowedTypes = onlyOfThisType == null ? null : EnumSet.of(onlyOfThisType);

            int n = 0;
            for (VariantContext vc : tracker.getAllVariantContexts(allowedTypes, context.getLocation(), onlyContextsStartinAtCurrentPosition, takeFirstOnly) ) {
                n++;
                if ( printContexts ) out.printf("       %s%n", vc);
            }

            if ( n > 0 && printContexts ) {
                out.printf("%s => had %d variant context objects%n", context.getLocation(), n);
                out.printf("---------------------------------------------%n");
            }
            
            return n;
        }
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer point, Integer sum) {
        return point + sum;
    }
}