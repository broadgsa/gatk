/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.examples;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;

/**
 * A simple Bayesian genotyper, that outputs a text based call format. Intended to be used only as an
 * example in the GATK publication.
 *
 * @author aaron
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_TOY, extraDocs = {CommandLineGATK.class} )
public class GATKPaperGenotyper extends LocusWalker<Integer,Long> implements TreeReducible<Long> {

    public static final double HUMAN_SNP_HETEROZYGOSITY = 1e-3;

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
        double likelihoods[] = getReferencePolarizedPriors(ref.getBase(),
                HUMAN_SNP_HETEROZYGOSITY,
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

        Integer sortedList[] = sortPermutation(likelihoods);

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

    private static Integer[] sortPermutation(final double[] A) {
        class comparator implements Comparator<Integer> {
            public int compare(Integer a, Integer b) {
                if (A[a.intValue()] < A[b.intValue()]) {
                    return -1;
                }
                if (A[a.intValue()] == A[b.intValue()]) {
                    return 0;
                }
                if (A[a.intValue()] > A[b.intValue()]) {
                    return 1;
                }
                return 0;
            }
        }
        Integer[] permutation = new Integer[A.length];
        for (int i = 0; i < A.length; i++) {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }

    /**
     * Takes reference base, and three priors for hom-ref, het, hom-var, and fills in the priors vector
     * appropriately.
     *
     * Suppose A is the reference base, and we are given the probability of being hom-ref, het, and hom-var,
     * and that pTriSateGenotype is the true probability of observing reference A and a true genotype of B/C
     * then this sets the priors to:
     *
     * AA = hom-ref
     * AC = AG = AT = (het - pTriStateGenotype) / 3
     * CC = GG = TT = hom-var / 3
     * CG = CT = GT = pTriStateGenotype / 3
     *
     * So that we get:
     *
     * hom-ref + 3 * (het - pTriStateGenotype) / 3 + 3 * hom-var / 3 + 3 * pTriStateGenotype
     * hom-ref + het - pTriStateGenotype + hom-var + pTriStateGenotype
     * hom-ref + het + hom-var
     * = 1
     *
     * @param ref
     * @param heterozyosity
     * @param pRefError
     */
    public static double[] getReferencePolarizedPriors(byte ref, double heterozyosity, double pRefError ) {
        if ( ! MathUtils.isBounded(pRefError, 0.0, 0.01) ) {
            throw new RuntimeException(String.format("BUG: p Reference error is out of bounds (0.0 - 0.01) is allow range %f", pRefError));
        }

        double pTriStateGenotype = heterozyosity * pRefError;
//        if ( pTriStateGenotype >= heterozyosity ) {
//            throw new RuntimeException(String.format("p Tristate genotype %f is greater than the heterozygosity %f", pTriStateGenotype, heterozyosity));
//        }

        double pHomRef = heterozygosity2HomRefProbability(heterozyosity);
        double pHet    = heterozygosity2HetProbability(heterozyosity);
        double pHomVar = heterozygosity2HomVarProbability(heterozyosity);

        if (MathUtils.compareDoubles(pHomRef + pHet + pHomVar, 1.0) != 0) {
            throw new RuntimeException(String.format("BUG: Prior probabilities don't sum to one => %f, %f, %f", pHomRef, pHet, pHomVar));
        }

        double[] priors = new double[DiploidGenotype.values().length];

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double POfG;

            final double nOnRefHets = 3;
            final double nOffRefHets = 3;
            final double nHomVars = 3;

            if ( g.isHomRef(ref) )      { POfG = pHomRef; }
            else if ( g.isHomVar(ref) ) { POfG = pHomVar / nHomVars; }
            else if ( g.isHetRef(ref) ) { POfG = (pHet - pTriStateGenotype ) / nOnRefHets; }
            else                        { POfG = pTriStateGenotype / nOffRefHets; }

            priors[g.ordinal()] = Math.log10(POfG);
        }

        return priors;
    }

    /**
     *
     * @param h
     * @return
     */
    public static double heterozygosity2HomRefProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        double v = 1.0 - (3.0 * h / 2.0);
        if (MathUtils.isNegative(v)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return v;
    }

    public static double heterozygosity2HetProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return h;
    }

    public static double heterozygosity2HomVarProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return h / 2.0;
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


