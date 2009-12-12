package org.broadinstitute.sting.playground.gatk.walkers.papergenotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotypePriors;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;


/**
 * A simple Bayesian genotyper, that outputs a text based call format. Intended to be used only as an
 * example in the GATK publication.
 * @author aaron
 * @help.summary A simple, naive Bayesian genotyper that is used as an example locus walker in the GATK paper. THIS IS NOT TO BE USED FOR ANY ANALYSIS
 */
public class GATKPaperGenotyper extends LocusWalker<SimpleCall, Integer> implements TreeReducible<Integer> {

    // the possible diploid genotype strings
    private static enum GENOTYPE { AA, AC, AG, AT, CC, CG, CT, GG, GT, TT }

    // where to write the genotyping data to
    @Argument(fullName = "call_location", shortName = "cl", doc = "File to which calls should be written", required = true)
    public PrintStream outputStream;

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

        ReadBackedPileup pileup = context.getPileup();
        double likelihoods[] = DiploidGenotypePriors.getReferencePolarizedPriors(ref.getBase(),
                                                                                 DiploidGenotypePriors.HUMAN_HETEROZYGOSITY,
                                                                                 DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);
        byte bases[] = pileup.getBases();
        byte quals[] = pileup.getQuals();
        for (GENOTYPE genotype : GENOTYPE.values())
            for (int index = 0; index < bases.length; index++) {
                if (quals[index] > 0) {
                    double epsilon = Math.pow(10, quals[index] / -10.0);
                    byte pileupBase = bases[index];
                    for (char genotypeBase : genotype.toString().toCharArray())
                        if (genotypeBase == pileupBase)
                            likelihoods[genotype.ordinal()] += Math.log10(0.5 * (1 - epsilon) + epsilon / 3);
                        else
                            likelihoods[genotype.ordinal()] += Math.log10(epsilon / 3);
                }
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
    public Integer reduceInit() {
        return 0;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(SimpleCall value, Integer sum) {
        if (value != null) outputStream.println(value.toString());
        return sum + 1;
    }

    /**
     * A composite, 'reduce of reduces' function.
     *
     * @param lhs 'left-most' portion of data in the composite reduce.
     * @param rhs 'right-most' portion of data in the composite reduce.
     * @return The composite reduce type.
     */
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    /**
     * when we finish traversing, close the result list
     * @param result the final reduce result
     */
    public void onTraversalDone(Integer result) {
        out.println("Simple Genotyper genotyped " + result + "Loci.");
    }
}


