/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext;

import org.broad.tribble.TribbleException;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.EnumMap;

public class GenotypeLikelihoods {
    public final static int MAX_PL = Short.MAX_VALUE;

    //
    // There are two objects here because we are lazy in creating both representations
    // for this object: a vector of log10 Probs and the PL phred-scaled string.  Supports
    // having one set during initializating, and dynamic creation of the other, if needed
    //
    private double[] log10Likelihoods = null;
    private String likelihoodsAsString_PLs = null;

    public final static GenotypeLikelihoods fromPLField(String PLs) {
        return new GenotypeLikelihoods(PLs);
    }

    public final static GenotypeLikelihoods fromGLField(String GLs) {
        return new GenotypeLikelihoods(parseDeprecatedGLString(GLs));
    }

    public final static GenotypeLikelihoods fromLog10Likelihoods(double[] log10Likelihoods) {
        return new GenotypeLikelihoods(log10Likelihoods);
    }

    public final static GenotypeLikelihoods fromPLs(final int[] pls) {
        return new GenotypeLikelihoods(PLsToGLs(pls));
    }

    //
    // You must use the factory methods now
    //
    private GenotypeLikelihoods(String asString) {
        likelihoodsAsString_PLs = asString;
    }

    private GenotypeLikelihoods(double[] asVector) {
        log10Likelihoods = asVector;
    }

    /**
     * Returns the genotypes likelihoods in negative log10 vector format.  pr{AA} = x, this
     * vector returns math.log10(x) for each of the genotypes.  Can return null if the
     * genotype likelihoods are "missing".
     *
     * @return
     */
    public double[] getAsVector() {
        // assumes one of the likelihoods vector or the string isn't null
        if ( log10Likelihoods == null ) {
            // make sure we create the GL string if it doesn't already exist
            log10Likelihoods = parsePLsIntoLikelihoods(likelihoodsAsString_PLs);
        }

        return log10Likelihoods;
    }

    public int[] getAsPLs() {
        final double[] GLs = getAsVector();
        return GLs == null ? null : GLsToPLs(GLs);
    }

    public String toString() {
        return getAsString();
    }

    public String getAsString() {
        if ( likelihoodsAsString_PLs == null ) {
            // todo -- should we accept null log10Likelihoods and set PLs as MISSING?
            if ( log10Likelihoods == null )
                throw new TribbleException("BUG: Attempted to get likelihoods as strings and neither the vector nor the string is set!");
            likelihoodsAsString_PLs = convertLikelihoodsToPLString(log10Likelihoods);
        }

        return likelihoodsAsString_PLs;
    }

    //Return genotype likelihoods as an EnumMap with Genotypes as keys and likelihoods as values
    //Returns null in case of missing likelihoods
    public EnumMap<Genotype.Type,Double> getAsMap(boolean normalizeFromLog10){
        //Make sure that the log10likelihoods are set
        double[] likelihoods = normalizeFromLog10 ? MathUtils.normalizeFromLog10(getAsVector()) : getAsVector();
        if(likelihoods == null)
            return null;
        EnumMap<Genotype.Type,Double> likelihoodsMap = new EnumMap<Genotype.Type, Double>(Genotype.Type.class);
        likelihoodsMap.put(Genotype.Type.HOM_REF,likelihoods[Genotype.Type.HOM_REF.ordinal()-1]);
        likelihoodsMap.put(Genotype.Type.HET,likelihoods[Genotype.Type.HET.ordinal()-1]);
        likelihoodsMap.put(Genotype.Type.HOM_VAR, likelihoods[Genotype.Type.HOM_VAR.ordinal() - 1]);
        return likelihoodsMap;
    }

    //Return the neg log10 Genotype Quality (GQ) for the given genotype
    //Returns Double.NEGATIVE_INFINITY in case of missing genotype
    public double getLog10GQ(Genotype.Type genotype){
        return getQualFromLikelihoods(genotype.ordinal() - 1 /* NO_CALL IS FIRST */, getAsVector());
    }

    public static double getQualFromLikelihoods(int iOfChoosenGenotype, double[] likelihoods){
        if(likelihoods == null)
            return Double.NEGATIVE_INFINITY;

        double qual = Double.NEGATIVE_INFINITY;
        for (int i=0; i < likelihoods.length; i++) {
            if (i==iOfChoosenGenotype)
                continue;
            if (likelihoods[i] >= qual)
                qual = likelihoods[i];
        }

        // qual contains now max(likelihoods[k]) for all k != bestGTguess
        qual = likelihoods[iOfChoosenGenotype] - qual;

        if (qual < 0) {
            // QUAL can be negative if the chosen genotype is not the most likely one individually.
            // In this case, we compute the actual genotype probability and QUAL is the likelihood of it not being the chosen one
            double[] normalized = MathUtils.normalizeFromLog10(likelihoods);
            double chosenGenotype = normalized[iOfChoosenGenotype];
            return Math.log10(1.0 - chosenGenotype);
        } else {
            // invert the size, as this is the probability of making an error
            return -1 * qual;
        }
    }

    private final static double[] parsePLsIntoLikelihoods(String likelihoodsAsString_PLs) {
        if ( !likelihoodsAsString_PLs.equals(VCFConstants.MISSING_VALUE_v4) ) {
            String[] strings = likelihoodsAsString_PLs.split(",");
            double[] likelihoodsAsVector = new double[strings.length];
            try {
                for ( int i = 0; i < strings.length; i++ ) {
                    likelihoodsAsVector[i] = Integer.parseInt(strings[i]) / -10.0;
                }
            } catch (NumberFormatException e) {
                throw new UserException.MalformedVCF("The GL/PL tag contains non-integer values: " + likelihoodsAsString_PLs);
            }
            return likelihoodsAsVector;
        } else
            return null;
    }

    /**
     * Back-compatibility function to read old style GL formatted genotype likelihoods in VCF format
     * @param GLString
     * @return
     */
    private final static double[] parseDeprecatedGLString(String GLString) {
        if ( !GLString.equals(VCFConstants.MISSING_VALUE_v4) ) {
            String[] strings = GLString.split(",");
            double[] likelihoodsAsVector = new double[strings.length];
            for ( int i = 0; i < strings.length; i++ ) {
                likelihoodsAsVector[i] = Double.parseDouble(strings[i]);
            }
            return likelihoodsAsVector;
        }

        return null;
    }

    private final static String convertLikelihoodsToPLString(final double[] GLs) {
        if ( GLs == null )
            return VCFConstants.MISSING_VALUE_v4;

        final StringBuilder s = new StringBuilder();
        boolean first = true;
        for ( final int pl : GLsToPLs(GLs) ) {
            if ( ! first )
                s.append(",");
            else
                first = false;

            s.append(pl);
        }

        return s.toString();
    }

    private final static int[] GLsToPLs(final double[] GLs) {
        final int[] pls = new int[GLs.length];
        final double adjust = maxPL(GLs);

        for ( int i = 0; i < GLs.length; i++ ) {
            pls[i] = (int)Math.round(Math.min(-10 * (GLs[i] - adjust), MAX_PL));
        }

        return pls;
    }

    private final static double maxPL(final double[] GLs) {
        double adjust = Double.NEGATIVE_INFINITY;
        for ( double l : GLs ) adjust = Math.max(adjust, l);
        return adjust;
    }

    private final static double[] PLsToGLs(final int pls[]) {
        double[] likelihoodsAsVector = new double[pls.length];
        for ( int i = 0; i < pls.length; i++ ) {
            likelihoodsAsVector[i] = pls[i] / -10.0;
        }
        return likelihoodsAsVector;
    }

//    // -------------------------------------------------------------------------------------
//    //
//    // List interface functions
//    //
//    // -------------------------------------------------------------------------------------
//
//    private final void notImplemented() {
//        throw new ReviewedStingException("BUG: code not implemented");
//    }
//
//    @Override public int size() { return getAsVector().length; }
//    @Override public Double get(final int i) { return getAsVector()[i];}
//    @Override public Double set(final int i, final Double aDouble) { return getAsVector()[i] = aDouble; }
//    @Override public boolean isEmpty() { return false; }
//    @Override public Iterator<Double> iterator() { return Arrays.asList(ArrayUtils.toObject(getAsVector())).iterator(); }
//    @Override public Object[] toArray() { return ArrayUtils.toObject(getAsVector()); }
//
//    // none of these are implemented
//    @Override public boolean contains(final Object o) { notImplemented(); return false; }
//    @Override public <T> T[] toArray(final T[] ts) { notImplemented(); return null; }
//    @Override public boolean add(final Double aDouble) { notImplemented(); return false; }
//    @Override public boolean remove(final Object o) {notImplemented(); return false; }
//    @Override public boolean containsAll(final Collection<?> objects) { notImplemented(); return false; }
//    @Override public boolean addAll(final Collection<? extends Double> doubles) { notImplemented(); return false; }
//    @Override public boolean addAll(final int i, final Collection<? extends Double> doubles) { notImplemented(); return false; }
//    @Override public boolean removeAll(final Collection<?> objects) { notImplemented(); return false; }
//    @Override public boolean retainAll(final Collection<?> objects) { notImplemented(); return false; }
//    @Override public void clear() { notImplemented(); }
//    @Override public void add(final int i, final Double aDouble) { notImplemented(); }
//    @Override public Double remove(final int i) { notImplemented(); return null; }
//    @Override public int indexOf(final Object o) { notImplemented(); return -1; }
//    @Override public int lastIndexOf(final Object o) { notImplemented(); return 0; }
//    @Override public ListIterator<Double> listIterator() { notImplemented(); return null; }
//    @Override public ListIterator<Double> listIterator(final int i) { notImplemented(); return null; }
//    @Override public List<Double> subList(final int i, final int i1) { notImplemented(); return null; }

    // -------------------------------------------------------------------------------------
    //
    // Static conversion utilities, going from GL/PL index to allele index and vice versa.
    //
    // -------------------------------------------------------------------------------------

    /*
     * Class representing the 2 alleles (or rather their indexes into VariantContext.getAllele()) corresponding to a specific PL index.
     * Note that the reference allele is always index=0.
     */
    public static class GenotypeLikelihoodsAllelePair {
        public final int alleleIndex1, alleleIndex2;

        public GenotypeLikelihoodsAllelePair(final int alleleIndex1, final int alleleIndex2) {
            this.alleleIndex1 = alleleIndex1;
            this.alleleIndex2 = alleleIndex2;
        }
    }

    /**
     * The maximum number of alleles that we can represent as genotype likelihoods
     */
    public final static int MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED = 50;

    /*
     * a cache of the PL index to the 2 alleles it represents over all possible numbers of alternate alleles
     */
    private final static GenotypeLikelihoodsAllelePair[] PLIndexToAlleleIndex = calculatePLcache(MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);

    private static GenotypeLikelihoodsAllelePair[] calculatePLcache(final int altAlleles) {
        final int numLikelihoods = calculateNumLikelihoods(1+altAlleles, 2);
        final GenotypeLikelihoodsAllelePair[] cache = new GenotypeLikelihoodsAllelePair[numLikelihoods];

        // for all possible combinations of 2 alleles
        for ( int allele1 = 0; allele1 <= altAlleles; allele1++ ) {
            for ( int allele2 = allele1; allele2 <= altAlleles; allele2++ ) {
                cache[calculatePLindex(allele1, allele2)] = new GenotypeLikelihoodsAllelePair(allele1, allele2);
            }
        }

        // a bit of sanity checking
        for ( int i = 0; i < cache.length; i++ ) {
            if ( cache[i] == null )
                throw new ReviewedStingException("BUG: cache entry " + i + " is unexpected null");
        }
            
        return cache;
    }

    /**
    * Compute how many likelihood elements are associated with the given number of alleles
    * Equivalent to asking in how many ways N non-negative integers can add up to P is S(N,P)
    * where P = ploidy (number of chromosomes) and N = total # of alleles.
    * Each chromosome can be in one single state (0,...,N-1) and there are P of them.
    * Naive solution would be to store N*P likelihoods, but this is not necessary because we can't distinguish chromosome states, but rather
    * only total number of alt allele counts in all chromosomes.
    *
    * For example, S(3,2) = 6: For alleles A,B,C, on a diploid organism we have six possible genotypes:
    * AA,AB,BB,AC,BC,CC.
    * Another way of expressing is with vector (#of A alleles, # of B alleles, # of C alleles)
    * which is then, for ordering above, (2,0,0), (1,1,0), (0,2,0), (1,1,0), (0,1,1), (0,0,2)
    * In general, for P=2 (regular biallelic), then S(N,2) = N*(N+1)/2
    *
    * Recursive implementation:
    *   S(N,P) = sum_{k=0}^P S(N-1,P-k)
    *  because if we have N integers, we can condition 1 integer to be = k, and then N-1 integers have to sum to P-K
    * With initial conditions
    *   S(N,1) = N  (only way to have N integers add up to 1 is all-zeros except one element with a one. There are N of these vectors)
    *   S(1,P) = 1 (only way to have 1 integer add to P is with that integer P itself).
    *
    *   @param  numAlleles      Number of alleles (including ref)
    *   @param  ploidy          Ploidy, or number of chromosomes in set
    *   @return    Number of likelihood elements we need to hold.
    */
    public static int calculateNumLikelihoods(final int numAlleles, final int ploidy) {

        // fast, closed form solution for diploid samples (most common use case)
        if (ploidy==2)
            return numAlleles*(numAlleles+1)/2;

        if (numAlleles == 1)
            return 1;
        else if (ploidy == 1)
            return numAlleles;

        int acc =0;
        for (int k=0; k <= ploidy; k++ )
            acc += calculateNumLikelihoods(numAlleles-1, ploidy-k);

        return acc;

    }

    // As per the VCF spec: "the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.
    // In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
    // Assumes that allele1Index < allele2Index
    public static int calculatePLindex(final int allele1Index, final int allele2Index) {
        return (allele2Index * (allele2Index+1) / 2) + allele1Index;
    }

    /**
     * get the allele index pair for the given PL
     *
     * @param PLindex   the PL index
     * @return the allele index pair
     */
    public static GenotypeLikelihoodsAllelePair getAllelePair(final int PLindex) {
        // make sure that we've cached enough data
        if ( PLindex >= PLIndexToAlleleIndex.length )
            throw new ReviewedStingException("GATK limitation: cannot genotype more than " + MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED + " alleles");

        return PLIndexToAlleleIndex[PLindex];
    }

    // An index conversion from the deprecated PL ordering to the new VCF-based ordering for up to 3 alternate alleles
    protected static int[] PLindexConversion = new int[]{0, 1, 3, 6, 2, 4, 7, 5, 8, 9};

    /**
     * get the allele index pair for the given PL using the deprecated PL ordering:
     *    AA,AB,AC,AD,BB,BC,BD,CC,CD,DD instead of AA,AB,BB,AC,BC,CC,AD,BD,CD,DD.
     * Although it's painful to keep this conversion around, our DiploidSNPGenotypeLikelihoods class uses the deprecated
     *    ordering and I know with certainty that external users have built code on top of it; changing it now would
     *    cause a whole lot of heartache for our collaborators, so for now at least there's a standard conversion method.
     * This method assumes at most 3 alternate alleles.
     *
     * @param PLindex   the PL index
     * @return the allele index pair
     */
    @Deprecated
    public static GenotypeLikelihoodsAllelePair getAllelePairUsingDeprecatedOrdering(final int PLindex) {
        return getAllelePair(PLindexConversion[PLindex]);
    }

    /**
     * get the PL indexes (AA, AB, BB) for the given allele pair; assumes allele1Index <= allele2Index.
     *
     * @param allele1Index    the index in VariantContext.getAllele() of the first allele
     * @param allele2Index    the index in VariantContext.getAllele() of the second allele
     * @return the PL indexes
     */
    public static int[] getPLIndecesOfAlleles(final int allele1Index, final int allele2Index) {
        
        final int[] indexes = new int[3];
        indexes[0] = calculatePLindex(allele1Index, allele1Index);
        indexes[1] = calculatePLindex(allele1Index, allele2Index);
        indexes[2] = calculatePLindex(allele2Index, allele2Index);
        return indexes;
    }
}
