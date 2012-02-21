/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {

    // private final static boolean DEBUG = false;

    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    protected ExactAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
    }

    public List<Allele> getLog10PNonRef(final VariantContext vc,
                                        final double[][] log10AlleleFrequencyPriors,
                                        final AlleleFrequencyCalculationResult result) {

        GenotypesContext GLs = vc.getGenotypes();
        List<Allele> alleles = vc.getAlleles();

        // don't try to genotype too many alternate alleles
        if ( vc.getAlternateAlleles().size() > MAX_ALTERNATE_ALLELES_TO_GENOTYPE ) {
            logger.warn("this tool is currently set to genotype at most " + MAX_ALTERNATE_ALLELES_TO_GENOTYPE + " alternate alleles in a given context, but the context at " + vc.getChr() + ":" + vc.getStart() + " has " + (vc.getAlternateAlleles().size()) + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

            alleles = new ArrayList<Allele>(MAX_ALTERNATE_ALLELES_TO_GENOTYPE + 1);
            alleles.add(vc.getReference());
            alleles.addAll(chooseMostLikelyAlternateAlleles(vc, MAX_ALTERNATE_ALLELES_TO_GENOTYPE));
            GLs = UnifiedGenotyperEngine.subsetAlleles(vc, alleles, false);
        }

        //linearExact(GLs, log10AlleleFrequencyPriors[0], log10AlleleFrequencyLikelihoods, log10AlleleFrequencyPosteriors);
        linearExactMultiAllelic(GLs, alleles.size() - 1, log10AlleleFrequencyPriors, result, false);

        return alleles;
    }

    private static final class LikelihoodSum implements Comparable<LikelihoodSum> {
        public double sum = 0.0;
        public Allele allele;

        public LikelihoodSum(Allele allele) { this.allele = allele; }

        public int compareTo(LikelihoodSum other) {
            final double diff = sum - other.sum;
            return ( diff < 0.0 ) ? 1 : (diff > 0.0 ) ? -1 : 0;
        }
    }

    private static final int PL_INDEX_OF_HOM_REF = 0;
    private static final List<Allele> chooseMostLikelyAlternateAlleles(VariantContext vc, int numAllelesToChoose) {
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ )
            likelihoodSums[i] = new LikelihoodSum(vc.getAlternateAllele(i));

        // make sure that we've cached enough data
        if ( numOriginalAltAlleles > UnifiedGenotyperEngine.PLIndexToAlleleIndex.length - 1 )
            UnifiedGenotyperEngine.calculatePLcache(numOriginalAltAlleles);

        // based on the GLs, find the alternate alleles with the most probability; sum the GLs for the most likely genotype
        final ArrayList<double[]> GLs = getGLs(vc.getGenotypes());
        for ( final double[] likelihoods : GLs ) {
            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            if ( PLindexOfBestGL != PL_INDEX_OF_HOM_REF ) {
                int[] alleles = UnifiedGenotyperEngine.PLIndexToAlleleIndex[numOriginalAltAlleles][PLindexOfBestGL];
                if ( alleles[0] != 0 )
                    likelihoodSums[alleles[0]-1].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
                // don't double-count it
                if ( alleles[1] != 0 && alleles[1] != alleles[0] )
                    likelihoodSums[alleles[1]-1].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
            }
        }

        // sort them by probability mass and choose the best ones
        Collections.sort(Arrays.asList(likelihoodSums));
        final ArrayList<Allele> bestAlleles = new ArrayList<Allele>(numAllelesToChoose);
        for ( int i = 0; i < numAllelesToChoose; i++ )
            bestAlleles.add(likelihoodSums[i].allele);

        final ArrayList<Allele> orderedBestAlleles = new ArrayList<Allele>(numAllelesToChoose);
        for ( Allele allele : vc.getAlternateAlleles() ) {
            if ( bestAlleles.contains(allele) )
                orderedBestAlleles.add(allele);
        }
        
        return orderedBestAlleles;
    }
    
    private static final ArrayList<double[]> getGLs(GenotypesContext GLs) {
        ArrayList<double[]> genotypeLikelihoods = new ArrayList<double[]>(GLs.size());

        genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                double[] gls = sample.getLikelihoods().getAsVector();

                if ( MathUtils.sum(gls) < UnifiedGenotyperEngine.SUM_GL_THRESH_NOCALL )
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }

    // -------------------------------------------------------------------------------------
    //
    // Multi-allelic implementation.
    //
    // -------------------------------------------------------------------------------------

    private static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first

    // a wrapper around the int array so that we can make it hashable
    private static final class ExactACcounts {

        private final int[] counts;
        private int hashcode = -1;

        public ExactACcounts(final int[] counts) {
            this.counts = counts;
        }

        public int[] getCounts() {
            return counts;
        }

        @Override
        public boolean equals(Object obj) {
            return (obj instanceof ExactACcounts) ? Arrays.equals(counts, ((ExactACcounts)obj).counts) : false;
        }

        @Override
        public int hashCode() {
            if ( hashcode == -1 )
                hashcode = Arrays.hashCode(counts);
            return hashcode;
        }

        @Override
        public String toString() {
            StringBuffer sb = new StringBuffer();
            sb.append(counts[0]);
            for ( int i = 1; i < counts.length; i++ ) {
                sb.append("/");
                sb.append(counts[i]);
            }
            return sb.toString();
        }
    }

    // This class represents a column in the Exact AC calculation matrix
    private static final class ExactACset {

        // the counts of the various alternate alleles which this column represents
        final ExactACcounts ACcounts;

        // the column of the matrix
        final double[] log10Likelihoods;

        // mapping of column index for those columns upon which this one depends to the index into the PLs which is used as the transition to this column;
        // for example, in the biallelic case, the transition from k=0 to k=1 would be AB while the transition to k=2 would be BB.
        final HashMap<ExactACcounts, Integer> ACsetIndexToPLIndex = new HashMap<ExactACcounts, Integer>();

        // to minimize memory consumption, we know we can delete any sets in this list because no further sets will depend on them
        final ArrayList<ExactACcounts> dependentACsetsToDelete = new ArrayList<ExactACcounts>();


        public ExactACset(final int size, final ExactACcounts ACcounts) {
            this.ACcounts = ACcounts;
            log10Likelihoods = new double[size];
        }

        // sum of all the non-reference alleles
        public int getACsum() {
            int sum = 0;
            for ( int count : ACcounts.getCounts() )
                sum += count;
            return sum;
        }

        public boolean equals(Object obj) {
            return (obj instanceof ExactACset) ? ACcounts.equals(((ExactACset)obj).ACcounts) : false;
        }
    }

    public static void linearExactMultiAllelic(final GenotypesContext GLs,
                                               final int numAlternateAlleles,
                                               final double[][] log10AlleleFrequencyPriors,
                                               final AlleleFrequencyCalculationResult result,
                                               final boolean preserveData) {

        // make sure the PL cache has been initialized
        if ( UnifiedGenotyperEngine.PLIndexToAlleleIndex == null )
            UnifiedGenotyperEngine.calculatePLcache(5);

        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        // queue of AC conformations to process
        final LinkedList<ExactACset> ACqueue = new LinkedList<ExactACset>();

        // mapping of ExactACset indexes to the objects
        final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<ExactACcounts, ExactACset>(numChr+1);

        // add AC=0 to the queue
        int[] zeroCounts = new int[numAlternateAlleles];
        ExactACset zeroSet = new ExactACset(numSamples+1, new ExactACcounts(zeroCounts));
        ACqueue.add(zeroSet);
        indexesToACset.put(zeroSet.ACcounts, zeroSet);

	    // optimization: create the temporary storage for computing L(j,k) just once
	    final int maxPossibleDependencies = numAlternateAlleles + (numAlternateAlleles * (numAlternateAlleles + 1) / 2) + 1;
	    final double[][] tempLog10ConformationLikelihoods = new double[numSamples+1][maxPossibleDependencies];
	    for ( int i = 0; i < maxPossibleDependencies; i++ )
	        tempLog10ConformationLikelihoods[0][i] = Double.NEGATIVE_INFINITY;

        // keep processing while we have AC conformations that need to be calculated
        double maxLog10L = Double.NEGATIVE_INFINITY;
        while ( !ACqueue.isEmpty() ) {
            // compute log10Likelihoods
            final ExactACset set = ACqueue.remove();
            final double log10LofKs = calculateAlleleCountConformation(set, genotypeLikelihoods, maxLog10L, numChr, preserveData, ACqueue, indexesToACset, log10AlleleFrequencyPriors, result, tempLog10ConformationLikelihoods);

            // adjust max likelihood seen if needed
            maxLog10L = Math.max(maxLog10L, log10LofKs);
        }
    }

    private static final class DependentSet {
        public final int[] ACcounts;
        public final int PLindex;
        
        public DependentSet(final int[] ACcounts, final int PLindex) {
            this.ACcounts = ACcounts;
            this.PLindex = PLindex;
        }
    }

    private static double calculateAlleleCountConformation(final ExactACset set,
                                                           final ArrayList<double[]> genotypeLikelihoods,
                                                           final double maxLog10L,
                                                           final int numChr,
                                                           final boolean preserveData,
                                                           final LinkedList<ExactACset> ACqueue,
                                                           final HashMap<ExactACcounts, ExactACset> indexesToACset,
                                                           final double[][] log10AlleleFrequencyPriors,
                                                           final AlleleFrequencyCalculationResult result,
                                                           final double[][] tempLog10ConformationLikelihoods) {

        //if ( DEBUG )
        //    System.out.printf(" *** computing LofK for set=%s%n", set.ACcounts);

        // compute the log10Likelihoods
        computeLofK(set, genotypeLikelihoods, indexesToACset, log10AlleleFrequencyPriors, result, tempLog10ConformationLikelihoods);

        // clean up memory
        if ( !preserveData ) {
            for ( ExactACcounts index : set.dependentACsetsToDelete ) {
                indexesToACset.remove(index);
                //if ( DEBUG )
                //    System.out.printf(" *** removing used set=%s after seeing final dependent set=%s%n", index, set.ACcounts);
            }
        }

        final double log10LofK = set.log10Likelihoods[set.log10Likelihoods.length-1];

        // can we abort early because the log10Likelihoods are so small?
        if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
            //if ( DEBUG )
            //    System.out.printf(" *** breaking early set=%s log10L=%.2f maxLog10L=%.2f%n", set.ACcounts, log10LofK, maxLog10L);

            // no reason to keep this data around because nothing depends on it
            if ( !preserveData )
                indexesToACset.remove(set.ACcounts);

            return log10LofK;
        }

        // iterate over higher frequencies if possible
        final int ACwiggle = numChr - set.getACsum();
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;

        final int numAltAlleles = set.ACcounts.getCounts().length;

        // genotype likelihoods are a linear vector that can be thought of as a row-wise upper triangular matrix of log10Likelihoods.
        // so e.g. with 2 alt alleles the likelihoods are AA,AB,AC,BB,BC,CC and with 3 alt alleles they are AA,AB,AC,AD,BB,BC,BD,CC,CD,DD.

        // add conformations for the k+1 case
        int PLindex = 0;
        for ( int allele = 0; allele < numAltAlleles; allele++ ) {
            final int[] ACcountsClone = set.ACcounts.getCounts().clone();
            ACcountsClone[allele]++;
            updateACset(ACcountsClone, numChr, set, ++PLindex, ACqueue, indexesToACset);
        }

        // add conformations for the k+2 case if it makes sense; note that the 2 new alleles may be the same or different
        if ( ACwiggle > 1 ) {
            final ArrayList<DependentSet> differentAlleles = new ArrayList<DependentSet>(numAltAlleles * numAltAlleles);
            final ArrayList<DependentSet> sameAlleles = new ArrayList<DependentSet>(numAltAlleles);

            for ( int allele_i = 0; allele_i < numAltAlleles; allele_i++ ) {
                for ( int allele_j = allele_i; allele_j < numAltAlleles; allele_j++ ) {
                    final int[] ACcountsClone = set.ACcounts.getCounts().clone();
                    ACcountsClone[allele_i]++;
                    ACcountsClone[allele_j]++;

                    if ( allele_i == allele_j )
                        sameAlleles.add(new DependentSet(ACcountsClone, ++PLindex));
                    else
                        differentAlleles.add(new DependentSet(ACcountsClone, ++PLindex));
                }
            }

            // IMPORTANT: we must first add the cases where the 2 new alleles are different so that the queue maintains its ordering
            for ( DependentSet dependent : differentAlleles )
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset);
            for ( DependentSet dependent : sameAlleles )
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset);
        }

        // determine which is the last dependent set in the queue (not necessarily the last one added above) so we can know when it is safe to clean up this column
        if ( !preserveData ) {
            final ExactACset lastSet = determineLastDependentSetInQueue(set.ACcounts, ACqueue);
            if ( lastSet != null )
                lastSet.dependentACsetsToDelete.add(set.ACcounts);
        }

        return log10LofK;
    }

    // adds the ExactACset represented by the ACcounts to the ACqueue if not already there (creating it if needed) and
    // also adds it as a dependency to the given callingSetIndex.
    // returns the ExactACset if that set was not already in the queue and null otherwise.
    private static void updateACset(final int[] ACcounts,
                                    final int numChr,
                                    final ExactACset callingSet,
                                    final int PLsetIndex,
                                    final Queue<ExactACset> ACqueue,
                                    final HashMap<ExactACcounts, ExactACset> indexesToACset) {
        final ExactACcounts index = new ExactACcounts(ACcounts);
        if ( !indexesToACset.containsKey(index) ) {
            ExactACset set = new ExactACset(numChr/2 +1, index);
            indexesToACset.put(index, set);
            ACqueue.add(set);
        }

        // add the given dependency to the set
        //if ( DEBUG )
        //    System.out.println(" *** adding dependency from " + index + " to " + callingSet.ACcounts);
        final ExactACset set = indexesToACset.get(index);
        set.ACsetIndexToPLIndex.put(callingSet.ACcounts, PLsetIndex);
    }

    private static ExactACset determineLastDependentSetInQueue(final ExactACcounts callingSetIndex, final LinkedList<ExactACset> ACqueue) {
        Iterator<ExactACset> reverseIterator = ACqueue.descendingIterator();
        while ( reverseIterator.hasNext() ) {
            final ExactACset queued = reverseIterator.next();
            if ( queued.ACsetIndexToPLIndex.containsKey(callingSetIndex) )
                return queued;
        }

        // shouldn't get here
        throw new ReviewedStingException("Error: no sets in the queue currently hold " + callingSetIndex + " as a dependent!");
    }

    private static void computeLofK(final ExactACset set,
                                    final ArrayList<double[]> genotypeLikelihoods,
                                    final HashMap<ExactACcounts, ExactACset> indexesToACset,
                                    final double[][] log10AlleleFrequencyPriors,
                                    final AlleleFrequencyCalculationResult result,
                                    final double[][] tempLog10ConformationLikelihoods) {

        set.log10Likelihoods[0] = 0.0; // the zero case
        final int totalK = set.getACsum();

        // special case for k = 0 over all k
        if ( totalK == 0 ) {
            for ( int j = 1; j < set.log10Likelihoods.length; j++ )
                set.log10Likelihoods[j] = set.log10Likelihoods[j-1] + genotypeLikelihoods.get(j)[HOM_REF_INDEX];
        }
        // k > 0 for at least one k
        else {
	        // deal with the non-AA possible conformations
            int conformationIndex = 1;
            for ( Map.Entry<ExactACcounts, Integer> mapping : set.ACsetIndexToPLIndex.entrySet() ) {
		        //if ( DEBUG )
		        //    System.out.printf(" *** evaluating set=%s which depends on set=%s%n", set.ACcounts, mapping.getKey());

                ExactACset dependent = indexesToACset.get(mapping.getKey());

                for ( int j = 1; j < set.log10Likelihoods.length; j++ ) {

                    if ( totalK <= 2*j ) { // skip impossible conformations
                        final double[] gl = genotypeLikelihoods.get(j);
                        tempLog10ConformationLikelihoods[j][conformationIndex] =
                                determineCoefficient(mapping.getValue(), j, set.ACcounts.getCounts(), totalK) + dependent.log10Likelihoods[j-1] + gl[mapping.getValue()];
                    } else {
                        tempLog10ConformationLikelihoods[j][conformationIndex] = Double.NEGATIVE_INFINITY;
                    }
                }

                conformationIndex++;
            }

	        // finally, deal with the AA case (which depends on previous cells in this column) and then update the L(j,k) value
            final int numPaths = set.ACsetIndexToPLIndex.size() + 1;
            for ( int j = 1; j < set.log10Likelihoods.length; j++ ) {

                if ( totalK < 2*j-1 ) {
                    final double[] gl = genotypeLikelihoods.get(j);
                    tempLog10ConformationLikelihoods[j][0] = MathUtils.log10Cache[2*j-totalK] + MathUtils.log10Cache[2*j-totalK-1] + set.log10Likelihoods[j-1] + gl[HOM_REF_INDEX];
                } else {
                    tempLog10ConformationLikelihoods[j][0] = Double.NEGATIVE_INFINITY;
                }

                final double logDenominator = MathUtils.log10Cache[2*j] + MathUtils.log10Cache[2*j-1];
                final double log10Max = MathUtils.approximateLog10SumLog10(tempLog10ConformationLikelihoods[j], numPaths);
                set.log10Likelihoods[j] = log10Max - logDenominator;
            }
        }

        final double log10LofK = set.log10Likelihoods[set.log10Likelihoods.length-1];

        // determine the power of theta to use
        int nonRefAlleles = 0;
        for ( int i = 0; i < set.ACcounts.getCounts().length; i++ ) {
            if ( set.ACcounts.getCounts()[i] > 0 )
                nonRefAlleles++;
        }

        // for k=0, we don't want to put that value into the likelihoods/posteriors matrix, but instead want to set the value in the results object
        if ( nonRefAlleles == 0 ) {
            result.log10LikelihoodOfAFzero = log10LofK;
            result.log10PosteriorOfAFzero = log10LofK + log10AlleleFrequencyPriors[0][0];
        } else {
            // update the likelihoods/posteriors vectors which are collapsed views of each of the various ACs
            for ( int i = 0; i < set.ACcounts.getCounts().length; i++ ) {
                int AC = set.ACcounts.getCounts()[i];
                result.log10AlleleFrequencyLikelihoods[i][AC] = MathUtils.approximateLog10SumLog10(result.log10AlleleFrequencyLikelihoods[i][AC], log10LofK);

                final double prior = log10AlleleFrequencyPriors[nonRefAlleles-1][AC];
                result.log10AlleleFrequencyPosteriors[i][AC] = MathUtils.approximateLog10SumLog10(result.log10AlleleFrequencyPosteriors[i][AC], log10LofK + prior);
            }
        }
    }

    private static double determineCoefficient(int PLindex, final int j, final int[] ACcounts, final int totalK) {

        // the closed form representation generalized for multiple alleles is as follows:
        // AA: (2j - totalK) * (2j - totalK - 1)
        // AB: 2k_b * (2j - totalK)
        // AC: 2k_c * (2j - totalK)
        // BB: k_b * (k_b - 1)
        // BC: 2 * k_b * k_c
        // CC: k_c * (k_c - 1)

        final int numAltAlleles = ACcounts.length;

        // the AX het case
        if ( PLindex <= numAltAlleles )
            return MathUtils.log10Cache[2*ACcounts[PLindex-1]] + MathUtils.log10Cache[2*j-totalK];

        // find the 2 alternate alleles that are represented by this PL index
        int[] alleles = UnifiedGenotyperEngine.PLIndexToAlleleIndex[numAltAlleles][PLindex];

        final int k_i = ACcounts[alleles[0]-1];  // subtract one because ACcounts doesn't consider the reference allele

        // the hom var case (e.g. BB, CC, DD)
        final double coeff;
        if ( alleles[0] == alleles[1] ) {
            coeff = MathUtils.log10Cache[k_i] + MathUtils.log10Cache[k_i - 1];
        }
        // the het non-ref case (e.g. BC, BD, CD)
        else {
            final int k_j = ACcounts[alleles[1]-1];
            coeff = MathUtils.log10Cache[2] + MathUtils.log10Cache[k_i] + MathUtils.log10Cache[k_j];
        }

        return coeff;
    }


    // -------------------------------------------------------------------------------------
    //
    // Deprecated bi-allelic ~O(N) implementation.  Kept here for posterity.
    //
    // -------------------------------------------------------------------------------------

    /**
     * A simple data structure that holds the current, prev, and prev->prev likelihoods vectors
     * for the exact model calculation
     */
/*
    private final static class ExactACCache {
        double[] kMinus2, kMinus1, kMinus0;

        private final static double[] create(int n) {
            return new double[n];
        }

        public ExactACCache(int n) {
            kMinus2 = create(n);
            kMinus1 = create(n);
            kMinus0 = create(n);
        }

        final public void rotate() {
            double[] tmp = kMinus2;
            kMinus2 = kMinus1;
            kMinus1 = kMinus0;
            kMinus0 = tmp;
        }

        final public double[] getkMinus2() {
            return kMinus2;
        }

        final public double[] getkMinus1() {
            return kMinus1;
        }

        final public double[] getkMinus0() {
            return kMinus0;
        }
    }

    public int linearExact(GenotypesContext GLs,
                           double[] log10AlleleFrequencyPriors,
                           double[][] log10AlleleFrequencyLikelihoods,
                           double[][] log10AlleleFrequencyPosteriors) {
        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        final ExactACCache logY = new ExactACCache(numSamples+1);
        logY.getkMinus0()[0] = 0.0; // the zero case

        double maxLog10L = Double.NEGATIVE_INFINITY;
        boolean done = false;
        int lastK = -1;

        for (int k=0; k <= numChr && ! done; k++ ) {
            final double[] kMinus0 = logY.getkMinus0();

            if ( k == 0 ) { // special case for k = 0
                for ( int j=1; j <= numSamples; j++ ) {
                    kMinus0[j] = kMinus0[j-1] + genotypeLikelihoods.get(j)[0];
                }
            } else { // k > 0
                final double[] kMinus1 = logY.getkMinus1();
                final double[] kMinus2 = logY.getkMinus2();

                for ( int j=1; j <= numSamples; j++ ) {
                    final double[] gl = genotypeLikelihoods.get(j);
                    final double logDenominator = MathUtils.log10Cache[2*j] + MathUtils.log10Cache[2*j-1];

                    double aa = Double.NEGATIVE_INFINITY;
                    double ab = Double.NEGATIVE_INFINITY;
                    if (k < 2*j-1)
                        aa = MathUtils.log10Cache[2*j-k] + MathUtils.log10Cache[2*j-k-1] + kMinus0[j-1] + gl[0];

                    if (k < 2*j)
                        ab = MathUtils.log10Cache[2*k] + MathUtils.log10Cache[2*j-k]+ kMinus1[j-1] + gl[1];

                    double log10Max;
                    if (k > 1) {
                        final double bb = MathUtils.log10Cache[k] + MathUtils.log10Cache[k-1] + kMinus2[j-1] + gl[2];
                        log10Max = approximateLog10SumLog10(aa, ab, bb);
                    } else {
                        // we know we aren't considering the BB case, so we can use an optimized log10 function
                        log10Max = approximateLog10SumLog10(aa, ab);
                    }

                    // finally, update the L(j,k) value
                    kMinus0[j] = log10Max - logDenominator;
                }
            }

            // update the posteriors vector
            final double log10LofK = kMinus0[numSamples];
            log10AlleleFrequencyLikelihoods[0][k] = log10LofK;
            log10AlleleFrequencyPosteriors[0][k] = log10LofK + log10AlleleFrequencyPriors[k];

            // can we abort early?
            lastK = k;
            maxLog10L = Math.max(maxLog10L, log10LofK);
            if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
                //if ( DEBUG ) System.out.printf("  *** breaking early k=%d log10L=%.2f maxLog10L=%.2f%n", k, log10LofK, maxLog10L);
                done = true;
            }

            logY.rotate();
        }

        return lastK;
    }

    final static double approximateLog10SumLog10(double a, double b, double c) {
        return approximateLog10SumLog10(approximateLog10SumLog10(a, b), c);
    }
*/

}
