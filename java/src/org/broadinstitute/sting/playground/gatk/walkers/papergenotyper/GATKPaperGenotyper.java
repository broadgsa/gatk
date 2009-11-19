package org.broadinstitute.sting.playground.gatk.walkers.papergenotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotypePriors;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.File;


/**
 * @author aaron
 *         <p/>
 *         Class GATKPaperGenotyper
 *         <p/>
 *         A simple Bayesian genotyper, that output a text based call format. Intended to be used only as an
 *         example in the GATK publication.
 */
public class GATKPaperGenotyper extends LocusWalker<SimpleCall, SimpleCallList> implements TreeReducible<SimpleCallList> {

    // the possible diploid genotype strings
    private static enum GENOTYPE { AA, AC, AG, AT, CC, CG, CT, GG, GT, TT }

    // the epsilon value we're using to model our error rate
    private static double EPSILON = 1e-4;

    @Argument(fullName = "call_location", shortName = "cl", doc = "File to which calls should be written", required = true)
    private File LOCATION = new File("genotyping.output");

    /**
     * our map function, which takes the reads spanning this locus, any associated reference ordered data,
     * and the reference information.  We output a simple genotype call as the result of this function
     *
     * @param tracker the reference ordered data tracker
     * @param ref     the reference information
     * @param context the locus context, which contains all of the read information
     * @return a SimpleCall, which stores the genotype we're calling and the LOD score
     */
    public SimpleCall map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

        ReadBackedPileup pileup = new ReadBackedPileup(context.getLocation(), ref.getBase(), context.getReads(), context.getOffsets());
        double likelihoods[] = DiploidGenotypePriors.getReferencePolarizedPriors(ref.getBase(),
                DiploidGenotypePriors.HUMAN_HETEROZYGOSITY,
                DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);

        for (GENOTYPE genotype : GENOTYPE.values())
            for (byte pileupBase : pileup.getBases()) {
                for (char genotypeBase : genotype.toString().toCharArray())
                    if (genotypeBase == pileupBase)
                        likelihoods[genotype.ordinal()] += 1 / 2 * (1 - EPSILON) + EPSILON / 3;
                    else
                        likelihoods[genotype.ordinal()] += EPSILON / 3;
            }

        Integer sortedList[] = Utils.SortPermutation(likelihoods);

        // get our reference genotype
        String refGenotype = (String.valueOf(ref.getBase()) + String.valueOf(ref.getBase())).toUpperCase();

        // create call using the best genotype (GENOTYPE.values()[sortedList[9]].toString())
        // and calculate the LOD score from best - ref (likelihoods[sortedList[9]] - likelihoods[sortedList[8])
        return new SimpleCall(context.getLocation(),
                GENOTYPE.values()[sortedList[9]].toString(),
                likelihoods[sortedList[9]] - likelihoods[GENOTYPE.valueOf(refGenotype).ordinal()]);
    }

    /**
     * Provide an initial value for reduce computations. In this case we simply return an empty list
     *
     * @return Initial value of reduce.
     */
    public SimpleCallList reduceInit() {
        return new SimpleCallList(LOCATION);
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public SimpleCallList reduce(SimpleCall value, SimpleCallList sum) {
        if (value != null) sum.add(value);
        return sum;
    }

    /**
     * A composite, 'reduce of reduces' function.
     *
     * @param lhs 'left-most' portion of data in the composite reduce.
     * @param rhs 'right-most' portion of data in the composite reduce.
     * @return The composite reduce type.
     */
    public SimpleCallList treeReduce(SimpleCallList lhs, SimpleCallList rhs) {
        lhs.addAll(rhs);
        return lhs;
    }

    /**
     * when we finish traversing, close the result list
     * @param result the final reduce result 
     */
    public void onTraversalDone(SimpleCallList result) {
        result.close();
    }
}


