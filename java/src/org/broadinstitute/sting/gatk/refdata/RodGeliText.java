/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.geli.GeliGenotypeCall;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class RodGeliText extends BasicReferenceOrderedDatum implements VariationRod, VariantBackedByGenotype {
    public enum Genotype_Strings {
        AA, AC, AG, AT, CC, CG, CT, GG, GT, TT
    }

    public GenomeLoc loc;
    public char refBase = 'N';
    public int depth;
    public int maxMappingQuality;
    public String bestGenotype = "NN";
    public double lodBtr;
    public double lodBtnb;
    public double[] genotypeLikelihoods = new double[10];

    public RodGeliText(final String name) {
        super(name);
    }

    public String delimiterRegex() {
        return "\\s+";
    }

    public boolean parseLine(Object header, String[] parts) throws IOException {
        if (parts.length < 18)
            throw new IOException("Invalid rodVariant row found -- too few elements.  Expected 18+, got " + parts.length);
        if (!parts[0].startsWith("#")) {
            loc = GenomeLocParser.createGenomeLoc(parts[0], Long.valueOf(parts[1]));
            refBase = Character.toUpperCase(parts[2].charAt(0));
            depth = Integer.valueOf(parts[3]);
            maxMappingQuality = Integer.valueOf(parts[4]);

            // UPPER case and sort
            char[] x = parts[5].toUpperCase().toCharArray();
            Arrays.sort(x);
            bestGenotype = new String(x);

            lodBtr = Double.valueOf(parts[6]);
            lodBtnb = Double.valueOf(parts[7]);

            for (int pieceIndex = 8, offset = 0; pieceIndex < 18; pieceIndex++, offset++) {
                genotypeLikelihoods[offset] = Double.valueOf(parts[pieceIndex]);
            }

            return true;
        }

        return false;
    }

    public String toString() {
        return String.format("%s\t%d\t%c\t%d\t%d\t%s\t%4.4f\t%4.4f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
                             loc.getContig(),
                             loc.getStart(),
                             refBase,
                             depth,
                             maxMappingQuality,
                             bestGenotype,
                             lodBtr,
                             lodBtnb,
                             genotypeLikelihoods[0],
                             genotypeLikelihoods[1],
                             genotypeLikelihoods[2],
                             genotypeLikelihoods[3],
                             genotypeLikelihoods[4],
                             genotypeLikelihoods[5],
                             genotypeLikelihoods[6],
                             genotypeLikelihoods[7],
                             genotypeLikelihoods[8],
                             genotypeLikelihoods[9]
        );
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    @Override
    public String getReference() {
        return String.valueOf(this.refBase);
    }


    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    @Override
    public double getNegLog10PError() {
        return Math.abs(lodBtr);
    }

    /**
     * gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference base.  This is returned as a string list with no guarantee ordering
     * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
     * frequency).
     *
     * @return an alternate allele list
     */
    @Override
    public List<String> getAlternateAlleleList() {
        List<String> list = new ArrayList<String>();
        for (char base : bestGenotype.toCharArray())
            if (base != refBase)
                list.add(String.valueOf(base));
        return list;
    }

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles.
     *
     * @return an alternate allele list
     */
    @Override
    public List<String> getAlleleList() {
        List<String> list = new ArrayList<String>();
        if (this.bestGenotype.contains(getReference())) list.add(getReference());
        for (char c : this.bestGenotype.toCharArray())
            if (c != Utils.stringToChar(getReference()))
                list.add(String.valueOf(c));
        return list;
    }

    public String getRefBasesFWD() {
        return String.format("%c", getRefSnpFWD());
    }

    public char getRefSnpFWD() throws IllegalStateException {
        return refBase;
    }

    public String getAltBasesFWD() {
        return String.format("%c", getAltSnpFWD());
    }

    public char getAltSnpFWD() throws IllegalStateException {
        // both ref and bestGenotype have been uppercased, so it's safe to use ==
        char c = (bestGenotype.charAt(0) == refBase) ? bestGenotype.charAt(1) : bestGenotype.charAt(0);
        //System.out.printf("%s : %c and %c%n", bestGenotype, refBase, c);
        return c;
    }

    public boolean isReference() {
        return refBase == bestGenotype.charAt(0) && refBase == bestGenotype.charAt(1);
    }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    @Override
    public double getNonRefAlleleFrequency() {
        return 1.0;
    }

    /** @return the VARIANT_TYPE of the current variant */
    @Override
    public VARIANT_TYPE getType() {
        return VARIANT_TYPE.SNP;
    }

    public boolean isSNP() {
        if (this.getReference().length() == 1)
            return (this.refBase != this.bestGenotype.charAt(0) || this.refBase != this.bestGenotype.charAt(1));
        return false;
    }

    public boolean isInsertion() {
        return false;
    }

    public boolean isDeletion() {
        return false;
    }

    public boolean isIndel() {
        return false;
    }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException in the case
     * of
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getAlternativeBaseForSNP() {
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        // we know that if we're a SNP, the alt is a single base
        if (this.bestGenotype.toString().charAt(0) == getReference().charAt(0))
            return this.bestGenotype.toString().charAt(1);
        return this.bestGenotype.toString().charAt(0);

    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getReferenceForSNP() {
        if (!isSNP()) throw new IllegalStateException("This site is not a SNP");
        // we know that if we're a SNP, the reference is a single base
        if (bestGenotype.toString().charAt(0) != getReference().charAt(0))
            return bestGenotype.toString().charAt(1);
        else
            return bestGenotype.toString().charAt(0);
    }

    public double getMAF() {
        return 0;
    }

    public double getHeterozygosity() {
        return 0;
    }

    public boolean isGenotype() {
        return true;
    }

    public double getVariationConfidence() {
        return lodBtr;
    }

    public double getConsensusConfidence() {
        return lodBtnb;
    }

    public List<String> getGenotype() throws IllegalStateException {
        return Arrays.asList(getBestGenotype());
    }

    public int getPloidy() throws IllegalStateException {
        return 2;
    }

    public boolean isBiallelic() {
        return true;
    }


    public int length() {
        return 1;
    }

    public char getReferenceBase() {
        return refBase;
    }

    public int getPileupDepth() {
        return depth;
    }

    public int getMaxMappingQuality() {
        return maxMappingQuality;
    }

    public String getBestGenotype() {
        return bestGenotype;
    }

    public double getLodBtr() {
        return lodBtr;
    }

    public double getLodBtnb() {
        return lodBtnb;
    }

    public double[] getGenotypeLikelihoods() {
        return genotypeLikelihoods;
    }

    public void adjustLikelihoods(double[] likelihoods) {
        for (int likelihoodIndex = 0; likelihoodIndex < likelihoods.length; likelihoodIndex++) {
            genotypeLikelihoods[likelihoodIndex] += likelihoods[likelihoodIndex];
        }

        String bestGenotype = "NN";
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        double nextBestLikelihood = Double.NEGATIVE_INFINITY;
        double refLikelihood = Double.NEGATIVE_INFINITY;

        for (int likelihoodIndex = 0; likelihoodIndex < likelihoods.length; likelihoodIndex++) {
            if (genotypeLikelihoods[likelihoodIndex] > bestLikelihood) {
                bestLikelihood = genotypeLikelihoods[likelihoodIndex];

                bestGenotype = Genotype_Strings.values()[likelihoodIndex].toString();
            }
        }

        for (int likelihoodIndex = 0; likelihoodIndex < likelihoods.length; likelihoodIndex++) {
            if (genotypeLikelihoods[likelihoodIndex] > nextBestLikelihood && genotypeLikelihoods[likelihoodIndex] < bestLikelihood) {
                nextBestLikelihood = genotypeLikelihoods[likelihoodIndex];
            }
        }

        for (int likelihoodIndex = 0; likelihoodIndex < likelihoods.length; likelihoodIndex++) {
            if (refBase == Genotype_Strings.values()[likelihoodIndex].toString().charAt(0) &&
                    refBase == Genotype_Strings.values()[likelihoodIndex].toString().charAt(1)) {
                refLikelihood = genotypeLikelihoods[likelihoodIndex];
            }
        }

        this.bestGenotype = bestGenotype;
        this.lodBtr = (bestLikelihood - refLikelihood);
        this.lodBtnb = (bestLikelihood - nextBestLikelihood);
    }


    /**
     * get the genotype
     *
     * @return a map in lexigraphical order of the genotypes
     */
    @Override
    public Genotype getCalledGenotype() {
        return new GeliGenotypeCall(refBase, getLocation(), bestGenotype, lodBtnb);
    }

    /**
     * get the likelihoods
     *
     * @return an array in lexigraphical order of the likelihoods
     */
    @Override
    public List<Genotype> getGenotypes() {
        List<Genotype> ret = new ArrayList<Genotype>();
        ret.add(new GeliGenotypeCall(refBase, getLocation(), bestGenotype, lodBtnb));
        return ret;
    }

    /**
     * do we have the specified genotype?  not all backedByGenotypes
     * have all the genotype data.
     *
     * @param x the genotype
     *
     * @return true if available, false otherwise
     */
    @Override
    public boolean hasGenotype(DiploidGenotype x) {
        return (x.toString().equals(this.getAltBasesFWD()));
    }
}
