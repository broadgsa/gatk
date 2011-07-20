/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.examples;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypePriors;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;


/**
 * A simple Bayesian genotyper, that outputs a text based call format. Intended to be used only as an
 * example in the GATK publication.
 *
 * @author aaron
 */
public class GATKPaperGenotyper extends LocusWalker<Integer,Long> implements TreeReducible<Long> {
    // the possible diploid genotype strings
    private static enum GENOTYPE { AA, AC, AG, AT, CC, CG, CT, GG, GT, TT }

    @Output
    private PrintStream out;

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
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

        ReadBackedPileup pileup = context.getBasePileup().getPileupWithoutMappingQualityZeroReads();
        double likelihoods[] = DiploidSNPGenotypePriors.getReferencePolarizedPriors(ref.getBase(),
                DiploidSNPGenotypePriors.HUMAN_HETEROZYGOSITY,
                0.01);
        // get the bases and qualities from the pileup
        byte bases[] = pileup.getBases();
        byte quals[] = pileup.getQuals();

        // for each genotype, determine it's likelihood value
        for (GENOTYPE genotype : GENOTYPE.values())
            for (int index = 0; index < bases.length; index++) {
                // our epsilon is the de-Phred scored base quality
                double epsilon = Math.pow(10, quals[index] / -10.0);

                byte pileupBase = bases[index];
                double p = 0;
                for (char r : genotype.toString().toCharArray())
                    p += r == pileupBase ? 1 - epsilon : epsilon / 3;
                likelihoods[genotype.ordinal()] += Math.log10(p / genotype.toString().length());
            }

        Integer sortedList[] = MathUtils.sortPermutation(likelihoods);

        // create call using the best genotype (GENOTYPE.values()[sortedList[9]].toString())
        // and calculate the LOD score from best - next best (9 and 8 in the sorted list, since the best likelihoods are closest to zero)
        GENOTYPE selectedGenotype = GENOTYPE.values()[sortedList[sortedList.length-1]];
        double lod = likelihoods[sortedList[sortedList.length-1]] - likelihoods[sortedList[sortedList.length-2]];

        if (lod > LODScore) {
            out.printf("%s\t%s\t%.4f\t%c%n", context.getLocation(), selectedGenotype, lod, (char)ref.getBase());
            return 1;
        }

        return 0;
    }

    /**
     * Provide an initial value for reduce computations. In this case we simply return an empty list
     *
     * @return Initial value of reduce.
     */
    public Long reduceInit() {
        return 0L;
    }

    /**
     * Outputs the number of genotypes called.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }

    /**
     * A composite, 'reduce of reduces' function.
     *
     * @param lhs 'left-most' portion of data in the composite reduce.
     * @param rhs 'right-most' portion of data in the composite reduce.
     * @return The composite reduce type.
     */
    public Long treeReduce(Long lhs, Long rhs) {
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


