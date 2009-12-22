package org.broadinstitute.sting.playground.gatk.walkers.diagnostics;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodGenotypeChipAsGFF;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.BaseUtils;

/**
 * Takes a BAM file and a Hapmap-chip file (via the -hc argument) and creates a table of reference allele
 * percentage and alternate allele percentage for het, homvar, and other genotypes.
 */
public class AlleleBalanceInspector extends LocusWalker<Integer, Integer> {
    private int item = 1;
    public void initialize() {
        out.printf("item\tlocus\tref\tgenotype\tstate\tdepth\trefdepth\taltdepth\trefpct\taltpct%n");
    }

    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        RodGenotypeChipAsGFF hc = (RodGenotypeChipAsGFF) tracker.lookup("hapmap-chip", null);

        return hc != null && hc.getCalledGenotype().isVariant(ref.getBase());
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        RodGenotypeChipAsGFF hc = (RodGenotypeChipAsGFF) tracker.lookup("hapmap-chip", null);

        String state;
        if (hc.getCalledGenotype().isHet()) {
            state = "het";
        } else if (hc.getCalledGenotype().isHom()) {
            state = "homvar";
        } else {
            state = "other";
        }

        int refIndex = ref.getBaseIndex();
        int altIndex = -1;
        for (char base : hc.getCalledGenotype().getBases().toCharArray()) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);

            if (baseIndex != refIndex) {
                altIndex = baseIndex;
            }
        }

        int[] baseCounts = context.getPileup().getBaseCounts();
        double sum = (double) (baseCounts[refIndex] + baseCounts[altIndex]);
        double refPct = ((double) baseCounts[refIndex])/sum;
        double altPct = ((double) baseCounts[altIndex])/sum;

        out.printf("%d\t%s\t%c\t%s\t%s\t%d\t%d\t%d\t%f\t%f%n",
                   item++,
                   context.getLocation(),
                   ref.getBase(),
                   hc.getCalledGenotype().getBases(),
                   state,
                   context.getPileup().getReads().size(),
                   baseCounts[refIndex],
                   baseCounts[altIndex], refPct, altPct);

        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
