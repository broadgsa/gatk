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

import java.util.EnumMap;

public class GenotypeLikelihoods {
    public static final boolean CAP_PLS = false;
    public static final int PL_CAP = 255;

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

    //
    // You must use the factory methods now
    //
    protected GenotypeLikelihoods(String asString) {
        likelihoodsAsString_PLs = asString;
    }

    protected GenotypeLikelihoods(double[] asVector) {
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
            for ( int i = 0; i < strings.length; i++ ) {
                likelihoodsAsVector[i] = Integer.parseInt(strings[i]) / -10.0;
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

    private final static String convertLikelihoodsToPLString(double[] GLs) {
        if ( GLs == null )
            return VCFConstants.MISSING_VALUE_v4;

        StringBuilder s = new StringBuilder();

        double adjust = Double.NEGATIVE_INFINITY;
        for ( double l : GLs ) adjust = Math.max(adjust, l);

        boolean first = true;
        for ( double l : GLs ) {
            if ( ! first )
                s.append(",");
            else
                first = false;

            long PL = Math.round(-10 * (l - adjust));
            if ( CAP_PLS )
                PL = Math.min(PL, PL_CAP);
            s.append(Long.toString(PL));
        }

        return s.toString();
    }


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

    /*
     * a cache of the PL index to the 2 alleles it represents over all possible numbers of alternate alleles
     */
    private static GenotypeLikelihoodsAllelePair[] PLIndexToAlleleIndex = new GenotypeLikelihoodsAllelePair[]{ new GenotypeLikelihoodsAllelePair(0, 0) };

    private static void calculatePLcache(final int minIndex) {
        // how many alternate alleles do we need to calculate for?
        int altAlleles = 0;
        int numLikelihoods = 1;
        while ( numLikelihoods <= minIndex ) {
            altAlleles++;
            numLikelihoods += altAlleles + 1;
        }

        PLIndexToAlleleIndex = new GenotypeLikelihoodsAllelePair[numLikelihoods];

        // for all possible combinations of 2 alleles
        for ( int allele1 = 0; allele1 <= altAlleles; allele1++ ) {
            for ( int allele2 = allele1; allele2 <= altAlleles; allele2++ ) {
                PLIndexToAlleleIndex[calculatePLindex(allele1, allele2)] = new GenotypeLikelihoodsAllelePair(allele1, allele2);
            }
        }
    }

    // how many likelihoods are associated with the given number of alternate alleles?
    public static int calculateNumLikelihoods(int numAltAlleles) {
        int numLikelihoods = 1;
        for ( int i = 1; i <= numAltAlleles; i++ )
            numLikelihoods += i + 1;
        return numLikelihoods;
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
            calculatePLcache(PLindex);

        return PLIndexToAlleleIndex[PLindex];
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
