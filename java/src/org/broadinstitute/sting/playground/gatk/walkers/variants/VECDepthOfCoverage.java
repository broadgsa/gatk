package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;

/**
 * Created by IntelliJ IDEA.
 * User: michaelmelgar
 * Date: Jun 22, 2009
 * Time: 6:04:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class VECDepthOfCoverage implements VariantExclusionCriterion {
    private int maximum = 200;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            maximum = Integer.valueOf(arguments);
        }
    }

    public boolean exclude(char ref, LocusContext context, rodVariants variant) {
        return context.getReads().size() > maximum;
    }
}
