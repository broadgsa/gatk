/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.variant;

import com.google.java.contract.Requires;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.variant.variantcontext.*;

import java.util.*;

public class GATKVariantContextUtils {

    public static final int DEFAULT_PLOIDY = 2;
    public static final double SUM_GL_THRESH_NOCALL = -0.1; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.
    private static final List<Allele> NO_CALL_ALLELES = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * create a genome location, given a variant context
     * @param genomeLocParser parser
     * @param vc the variant context
     * @return the genomeLoc
     */
    public static final GenomeLoc getLocation(GenomeLocParser genomeLocParser,VariantContext vc) {
        return genomeLocParser.createGenomeLoc(vc.getChr(), vc.getStart(), vc.getEnd(), true);
    }

    /**
     * Returns true iff VC is an non-complex indel where every allele represents an expansion or
     * contraction of a series of identical bases in the reference.
     *
     * For example, suppose the ref bases are CTCTCTGA, which includes a 3x repeat of CTCTCT
     *
     * If VC = -/CT, then this function returns true because the CT insertion matches exactly the
     * upcoming reference.
     * If VC = -/CTA then this function returns false because the CTA isn't a perfect match
     *
     * Now consider deletions:
     *
     * If VC = CT/- then again the same logic applies and this returns true
     * The case of CTA/- makes no sense because it doesn't actually match the reference bases.
     *
     * The logic of this function is pretty simple.  Take all of the non-null alleles in VC.  For
     * each insertion allele of n bases, check if that allele matches the next n reference bases.
     * For each deletion allele of n bases, check if this matches the reference bases at n - 2 n,
     * as it must necessarily match the first n bases.  If this test returns true for all
     * alleles you are a tandem repeat, otherwise you are not.
     *
     * @param vc
     * @param refBasesStartingAtVCWithPad not this is assumed to include the PADDED reference
     * @return
     */
    @Requires({"vc != null", "refBasesStartingAtVCWithPad != null && refBasesStartingAtVCWithPad.length > 0"})
    public static boolean isTandemRepeat(final VariantContext vc, final byte[] refBasesStartingAtVCWithPad) {
        final String refBasesStartingAtVCWithoutPad = new String(refBasesStartingAtVCWithPad).substring(1);
        if ( ! vc.isIndel() ) // only indels are tandem repeats
            return false;

        final Allele ref = vc.getReference();

        for ( final Allele allele : vc.getAlternateAlleles() ) {
            if ( ! isRepeatAllele(ref, allele, refBasesStartingAtVCWithoutPad) )
                return false;
        }

        // we've passed all of the tests, so we are a repeat
        return true;
    }

    /**
     *
     * @param vc
     * @param refBasesStartingAtVCWithPad
     * @return
     */
    @Requires({"vc != null", "refBasesStartingAtVCWithPad != null && refBasesStartingAtVCWithPad.length > 0"})
    public static Pair<List<Integer>,byte[]> getNumTandemRepeatUnits(final VariantContext vc, final byte[] refBasesStartingAtVCWithPad) {
        final boolean VERBOSE = false;
        final String refBasesStartingAtVCWithoutPad = new String(refBasesStartingAtVCWithPad).substring(1);
        if ( ! vc.isIndel() ) // only indels are tandem repeats
            return null;

        final Allele refAllele = vc.getReference();
        final byte[] refAlleleBases = Arrays.copyOfRange(refAllele.getBases(), 1, refAllele.length());

        byte[] repeatUnit = null;
        final ArrayList<Integer> lengths = new ArrayList<Integer>();

        for ( final Allele allele : vc.getAlternateAlleles() ) {
            Pair<int[],byte[]> result = getNumTandemRepeatUnits(refAlleleBases, Arrays.copyOfRange(allele.getBases(), 1, allele.length()), refBasesStartingAtVCWithoutPad.getBytes());

            final int[] repetitionCount = result.first;
            // repetition count = 0 means allele is not a tandem expansion of context
            if (repetitionCount[0] == 0 || repetitionCount[1] == 0)
                return null;

            if (lengths.size() == 0) {
                lengths.add(repetitionCount[0]); // add ref allele length only once
            }
            lengths.add(repetitionCount[1]);  // add this alt allele's length

            repeatUnit = result.second;
            if (VERBOSE) {
                System.out.println("RefContext:"+refBasesStartingAtVCWithoutPad);
                System.out.println("Ref:"+refAllele.toString()+" Count:" + String.valueOf(repetitionCount[0]));
                System.out.println("Allele:"+allele.toString()+" Count:" + String.valueOf(repetitionCount[1]));
                System.out.println("RU:"+new String(repeatUnit));
            }
        }

        return new Pair<List<Integer>, byte[]>(lengths,repeatUnit);
    }

    protected static Pair<int[],byte[]> getNumTandemRepeatUnits(final byte[] refBases, final byte[] altBases, final byte[] remainingRefContext) {
         /* we can't exactly apply same logic as in basesAreRepeated() to compute tandem unit and number of repeated units.
           Consider case where ref =ATATAT and we have an insertion of ATAT. Natural description is (AT)3 -> (AT)5.
         */

        byte[] longB;
        // find first repeat unit based on either ref or alt, whichever is longer
        if (altBases.length > refBases.length)
            longB = altBases;
        else
            longB = refBases;

        // see if non-null allele (either ref or alt, whichever is longer) can be decomposed into several identical tandem units
        // for example, -*,CACA needs to first be decomposed into (CA)2
        final int repeatUnitLength = findRepeatedSubstring(longB);
        final byte[] repeatUnit = Arrays.copyOf(longB, repeatUnitLength);

        final int[] repetitionCount = new int[2];
//        repetitionCount[0] = findNumberofRepetitions(repeatUnit, ArrayUtils.addAll(refBases, remainingRefContext));
//        repetitionCount[1] = findNumberofRepetitions(repeatUnit, ArrayUtils.addAll(altBases, remainingRefContext));
        int repetitionsInRef = findNumberofRepetitions(repeatUnit,refBases);
        repetitionCount[0] = findNumberofRepetitions(repeatUnit, ArrayUtils.addAll(refBases, remainingRefContext))-repetitionsInRef;
        repetitionCount[1] = findNumberofRepetitions(repeatUnit, ArrayUtils.addAll(altBases, remainingRefContext))-repetitionsInRef;

        return new Pair<int[], byte[]>(repetitionCount, repeatUnit);

    }

    /**
     * Find out if a string can be represented as a tandem number of substrings.
     * For example ACTACT is a 2-tandem of ACT,
     * but ACTACA is not.
     *
     * @param bases                 String to be tested
     * @return                      Length of repeat unit, if string can be represented as tandem of substring (if it can't
     *                              be represented as one, it will be just the length of the input string)
     */
    public static int findRepeatedSubstring(byte[] bases) {

        int repLength;
        for (repLength=1; repLength <=bases.length; repLength++) {
            final byte[] candidateRepeatUnit = Arrays.copyOf(bases,repLength);
            boolean allBasesMatch = true;
            for (int start = repLength; start < bases.length; start += repLength ) {
                // check that remaining of string is exactly equal to repeat unit
                final byte[] basePiece = Arrays.copyOfRange(bases,start,start+candidateRepeatUnit.length);
                if (!Arrays.equals(candidateRepeatUnit, basePiece)) {
                    allBasesMatch = false;
                    break;
                }
            }
            if (allBasesMatch)
                return repLength;
        }

        return repLength;
    }

    /**
     * Helper routine that finds number of repetitions a string consists of.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnit             Substring
     * @param testString             String to test
     * @return                       Number of repetitions (0 if testString is not a concatenation of n repeatUnit's
     */
    public static int findNumberofRepetitions(byte[] repeatUnit, byte[] testString) {
        int numRepeats = 0;
        for (int start = 0; start < testString.length; start += repeatUnit.length) {
            int end = start + repeatUnit.length;
            byte[] unit = Arrays.copyOfRange(testString,start, end);
            if(Arrays.equals(unit,repeatUnit))
                numRepeats++;
            else
                return numRepeats;
        }
        return numRepeats;
    }

    /**
     * Helper function for isTandemRepeat that checks that allele matches somewhere on the reference
     * @param ref
     * @param alt
     * @param refBasesStartingAtVCWithoutPad
     * @return
     */
    protected static boolean isRepeatAllele(final Allele ref, final Allele alt, final String refBasesStartingAtVCWithoutPad) {
        if ( ! Allele.oneIsPrefixOfOther(ref, alt) )
            return false; // we require one allele be a prefix of another

        if ( ref.length() > alt.length() ) { // we are a deletion
            return basesAreRepeated(ref.getBaseString(), alt.getBaseString(), refBasesStartingAtVCWithoutPad, 2);
        } else { // we are an insertion
            return basesAreRepeated(alt.getBaseString(), ref.getBaseString(), refBasesStartingAtVCWithoutPad, 1);
        }
    }

    protected static boolean basesAreRepeated(final String l, final String s, final String ref, final int minNumberOfMatches) {
        final String potentialRepeat = l.substring(s.length()); // skip s bases

        for ( int i = 0; i < minNumberOfMatches; i++) {
            final int start = i * potentialRepeat.length();
            final int end = (i+1) * potentialRepeat.length();
            if ( ref.length() < end )
                 return false; // we ran out of bases to test
            final String refSub = ref.substring(start, end);
            if ( ! refSub.equals(potentialRepeat) )
                return false; // repeat didn't match, fail
        }

        return true; // we passed all tests, we matched
    }

    /**
     * Assign genotypes (GTs) to the samples in the Variant Context greedily based on the PLs
     *
     * @param vc            variant context with genotype likelihoods
     * @return genotypes context
     */
    public static GenotypesContext assignDiploidGenotypes(final VariantContext vc) {
        return subsetDiploidAlleles(vc, vc.getAlleles(), true);
    }

    /**
     * Split variant context into its biallelic components if there are more than 2 alleles
     *
     * For VC has A/B/C alleles, returns A/B and A/C contexts.
     * Genotypes are all no-calls now (it's not possible to fix them easily)
     * Alleles are right trimmed to satisfy VCF conventions
     *
     * If vc is biallelic or non-variant it is just returned
     *
     * Chromosome counts are updated (but they are by definition 0)
     *
     * @param vc a potentially multi-allelic variant context
     * @return a list of bi-allelic (or monomorphic) variant context
     */
    public static List<VariantContext> splitVariantContextToBiallelics(final VariantContext vc) {
        if ( ! vc.isVariant() || vc.isBiallelic() )
            // non variant or biallelics already satisfy the contract
            return Collections.singletonList(vc);
        else {
            final List<VariantContext> biallelics = new LinkedList<VariantContext>();

            for ( final Allele alt : vc.getAlternateAlleles() ) {
                VariantContextBuilder builder = new VariantContextBuilder(vc);
                final List<Allele> alleles = Arrays.asList(vc.getReference(), alt);
                builder.alleles(alleles);
                builder.genotypes(subsetDiploidAlleles(vc, alleles, false));
                VariantContextUtils.calculateChromosomeCounts(builder, true);
                biallelics.add(reverseTrimAlleles(builder.make()));
            }

            return biallelics;
        }
    }

    /**
     * subset the Variant Context to the specific set of alleles passed in (pruning the PLs appropriately)
     *
     * @param vc                 variant context with genotype likelihoods
     * @param allelesToUse       which alleles from the vc are okay to use; *** must be in the same relative order as those in the original VC ***
     * @param assignGenotypes    true if we should update the genotypes based on the (subsetted) PLs
     * @return genotypes
     */
    public static GenotypesContext subsetDiploidAlleles(final VariantContext vc,
                                                 final List<Allele> allelesToUse,
                                                 final boolean assignGenotypes) {

        // the genotypes with PLs
        final GenotypesContext oldGTs = vc.getGenotypes();

        // samples
        final List<String> sampleIndices = oldGTs.getSampleNamesOrderedByName();

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create();

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // which PLs should be carried forward?
        ArrayList<Integer> likelihoodIndexesToUse = null;

        // an optimization: if we are supposed to use all (or none in the case of a ref call) of the alleles,
        // then we can keep the PLs as is; otherwise, we determine which ones to keep
        if ( numNewAltAlleles != numOriginalAltAlleles && numNewAltAlleles > 0 ) {
            likelihoodIndexesToUse = new ArrayList<Integer>(30);

            final boolean[] altAlleleIndexToUse = new boolean[numOriginalAltAlleles];
            for ( int i = 0; i < numOriginalAltAlleles; i++ ) {
                if ( allelesToUse.contains(vc.getAlternateAllele(i)) )
                    altAlleleIndexToUse[i] = true;
            }

            // numLikelihoods takes total # of alleles. Use default # of chromosomes (ploidy) = 2
            final int numLikelihoods = GenotypeLikelihoods.numLikelihoods(1 + numOriginalAltAlleles, DEFAULT_PLOIDY);
            for ( int PLindex = 0; PLindex < numLikelihoods; PLindex++ ) {
                final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindex);
                // consider this entry only if both of the alleles are good
                if ( (alleles.alleleIndex1 == 0 || altAlleleIndexToUse[alleles.alleleIndex1 - 1]) && (alleles.alleleIndex2 == 0 || altAlleleIndexToUse[alleles.alleleIndex2 - 1]) )
                    likelihoodIndexesToUse.add(PLindex);
            }
        }

        // create the new genotypes
        for ( int k = 0; k < oldGTs.size(); k++ ) {
            final Genotype g = oldGTs.get(sampleIndices.get(k));
            if ( !g.hasLikelihoods() ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), NO_CALL_ALLELES));
                continue;
            }

            // create the new likelihoods array from the alleles we are allowed to use
            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            double[] newLikelihoods;
            if ( likelihoodIndexesToUse == null ) {
                newLikelihoods = originalLikelihoods;
            } else {
                newLikelihoods = new double[likelihoodIndexesToUse.size()];
                int newIndex = 0;
                for ( int oldIndex : likelihoodIndexesToUse )
                    newLikelihoods[newIndex++] = originalLikelihoods[oldIndex];

                // might need to re-normalize
                newLikelihoods = MathUtils.normalizeFromLog10(newLikelihoods, false, true);
            }

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if ( MathUtils.sum(newLikelihoods) > SUM_GL_THRESH_NOCALL ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), NO_CALL_ALLELES));
            }
            else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);

                if ( numNewAltAlleles == 0 )
                    gb.noPL();
                else
                    gb.PL(newLikelihoods);

                // if we weren't asked to assign a genotype, then just no-call the sample
                if ( !assignGenotypes || MathUtils.sum(newLikelihoods) > SUM_GL_THRESH_NOCALL ) {
                    gb.alleles(NO_CALL_ALLELES);
                }
                else {
                    // find the genotype with maximum likelihoods
                    int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);
                    GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindex);

                    gb.alleles(Arrays.asList(allelesToUse.get(alleles.alleleIndex1), allelesToUse.get(alleles.alleleIndex2)));
                    if ( numNewAltAlleles != 0 ) gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
                }
                newGTs.add(gb.make());
            }
        }

        return newGTs;
    }

    public static VariantContext reverseTrimAlleles( final VariantContext inputVC ) {

        // see whether we need to trim common reference base from all alleles
        final int trimExtent = computeReverseClipping(inputVC.getAlleles(), inputVC.getReference().getDisplayString().getBytes(), 0, false);
        if ( trimExtent <= 0 || inputVC.getAlleles().size() <= 1 )
            return inputVC;

        final List<Allele> alleles = new ArrayList<Allele>();
        final GenotypesContext genotypes = GenotypesContext.create();
        final Map<Allele, Allele> originalToTrimmedAlleleMap = new HashMap<Allele, Allele>();

        for (final Allele a : inputVC.getAlleles()) {
            if (a.isSymbolic()) {
                alleles.add(a);
                originalToTrimmedAlleleMap.put(a, a);
            } else {
                // get bases for current allele and create a new one with trimmed bases
                final byte[] newBases = Arrays.copyOfRange(a.getBases(), 0, a.length()-trimExtent);
                final Allele trimmedAllele = Allele.create(newBases, a.isReference());
                alleles.add(trimmedAllele);
                originalToTrimmedAlleleMap.put(a, trimmedAllele);
            }
        }

        // now we can recreate new genotypes with trimmed alleles
        for ( final Genotype genotype : inputVC.getGenotypes() ) {
            final List<Allele> originalAlleles = genotype.getAlleles();
            final List<Allele> trimmedAlleles = new ArrayList<Allele>();
            for ( final Allele a : originalAlleles ) {
                if ( a.isCalled() )
                    trimmedAlleles.add(originalToTrimmedAlleleMap.get(a));
                else
                    trimmedAlleles.add(Allele.NO_CALL);
            }
            genotypes.add(new GenotypeBuilder(genotype).alleles(trimmedAlleles).make());
        }

        return new VariantContextBuilder(inputVC).stop(inputVC.getStart() + alleles.get(0).length() - 1).alleles(alleles).genotypes(genotypes).make();
    }

    public static int computeReverseClipping(final List<Allele> unclippedAlleles,
                                             final byte[] ref,
                                             final int forwardClipping,
                                             final boolean allowFullClip) {
        int clipping = 0;
        boolean stillClipping = true;

        while ( stillClipping ) {
            for ( final Allele a : unclippedAlleles ) {
                if ( a.isSymbolic() )
                    continue;

                // we need to ensure that we don't reverse clip out all of the bases from an allele because we then will have the wrong
                // position set for the VariantContext (although it's okay to forward clip it all out, because the position will be fine).
                if ( a.length() - clipping == 0 )
                    return clipping - (allowFullClip ? 0 : 1);

                if ( a.length() - clipping <= forwardClipping || a.length() - forwardClipping == 0 ) {
                    stillClipping = false;
                }
                else if ( ref.length == clipping ) {
                    if ( allowFullClip )
                        stillClipping = false;
                    else
                        return -1;
                }
                else if ( a.getBases()[a.length()-clipping-1] != ref[ref.length-clipping-1] ) {
                    stillClipping = false;
                }
            }
            if ( stillClipping )
                clipping++;
        }

        return clipping;
    }
}
