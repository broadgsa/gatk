package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;

import java.io.PrintStream;

/**
 * Walks over the input data set, calculating the percentage of loci covered by at least 20 reads.
 */
public class Percent20xCoverage extends LocusWalker<Integer, Long> implements TreeReducible<Long> {
    @Output(doc="Write count to this file instead of STDOUT")
    PrintStream out;

    @Argument(fullName="target_coverage", shortName="coverage", doc="Set the target coverage", required=false)
    private int targetCoverage = 20;


    private long totalLoci;
    private long totalCoverage;

    public void initialize () {
        totalLoci = 0;
        totalCoverage = 0;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        totalLoci++;
        int coverage = context.getBasePileup().size();
        totalCoverage += coverage;
        if (coverage >= targetCoverage)
            return 1;
        return 0;
    }

    public Long reduceInit() {
        return 0l;
    }

    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }

    public Long treeReduce(Long lhs, Long rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone( Long c ) {
        out.println(Math.floor(totalCoverage/totalLoci) + "\t" + (double) c/totalLoci);
    }
}
