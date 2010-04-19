/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers.papergenotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotypePriors;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;


/**
 * A simple Bayesian genotyper, that outputs a text based call format. Intended to be used only as an
 * example in the GATK publication.
 *
 * @author aaron
 */
public class GATKPaperGenotyper extends LocusWalker<SimpleCall, Integer> implements TreeReducible<Integer> {

    // the possible diploid genotype strings
    private static enum GENOTYPE { AA, AC, AG, AT, CC, CG, CT, GG, GT, TT }

    @Argument(fullName = "call_location", shortName = "cl", doc = "File to which calls should be written", required = true)
    private PrintStream outputStream;

    @Argument(fullName = "log_odds_score", shortName = "LOD", doc = "The LOD threshold for us to call confidently a genotype", required = false)
    private double LODScore = 3.0;

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
                                                                                 0.01);
        // get the bases and qualities from the pileup
        byte bases[] = pileup.getBases();
        byte quals[] = pileup.getQuals();

        // for each genotype, determine it's likelihood value
        for (GENOTYPE genotype : GENOTYPE.values())
            for (int index = 0; index < bases.length; index++) {
                if (quals[index] > 0) {
                    // our epsilon is the de-Phred scored base quality
                    double epsilon = Math.pow(10, quals[index] / -10.0);

                    byte pileupBase = bases[index];
                    double p = 0;
                    for (char r : genotype.toString().toCharArray())
                        p += r == pileupBase ? 1 - epsilon : epsilon / 3;
                    likelihoods[genotype.ordinal()] += Math.log10(p / genotype.toString().length());
                }
            }

        Integer sortedList[] = MathUtils.sortPermutation(likelihoods);

        // create call using the best genotype (GENOTYPE.values()[sortedList[9]].toString())
        // and calculate the LOD score from best - next best (9 and 8 in the sorted list, since the best likelihoods are closest to zero)
        return new SimpleCall(context.getLocation(),
                              GENOTYPE.values()[sortedList[9]].toString(),
                              likelihoods[sortedList[9]] - likelihoods[sortedList[8]],
                              ref.getBase());
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
     * Reduces a single map with the accumulator provided as the ReduceType. We filter out calls,
     * first making sure that the call is != null, secondly that the LOD score is above a moderate
     * threshold (in this case 3).
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(SimpleCall value, Integer sum) {
        if (value != null && value.LOD > LODScore) outputStream.println(value.toString());
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


