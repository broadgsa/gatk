package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.io.IOException;
import java.util.List;
import java.util.Arrays;

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

public class rodVariants extends BasicReferenceOrderedDatum implements AllelicVariant {
    private enum Genotype { AA, AC, AG, AT, CC, CG, CT, GG, GT, TT }

    private GenomeLoc loc;
    private char refBase = 'N';
    private int depth;
    private int maxMappingQuality;
    private String bestGenotype = "NN";
    private float lodBtr;
    private float lodBtnb;
    private float[] genotypeLikelihoods = new float[10];
    
    public rodVariants(final String name) { super(name); }

    public String delimiterRegex() { return "\\s+"; }

    public boolean parseLine(Object header, String[] parts) throws IOException {
        if (!parts[0].startsWith("#")) {
            loc = GenomeLocParser.createGenomeLoc(parts[0], Long.valueOf(parts[1]));
            refBase = parts[2].charAt(0);
            depth = Integer.valueOf(parts[3]);
            maxMappingQuality = Integer.valueOf(parts[4]);
            bestGenotype = parts[5];
            lodBtr = Float.valueOf(parts[6]);
            lodBtnb = Float.valueOf(parts[7]);

            for (int pieceIndex = 8, offset = 0; pieceIndex < 18; pieceIndex++, offset++) {
                genotypeLikelihoods[offset] = Float.valueOf(parts[pieceIndex]);
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

    public GenomeLoc getLocation() { return loc; }

    public char getReferenceBase() { return refBase; }

    public int getPileupDepth() { return depth; }

    public int getMaxMappingQuality() { return maxMappingQuality; }

    public String getBestGenotype() { return bestGenotype; }

    public float getLodBtr() { return lodBtr; }

    public float getLodBtnb() { return lodBtnb; }

    public float[] getGenotypeLikelihoods() { return genotypeLikelihoods; }

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

                bestGenotype = Genotype.values()[likelihoodIndex].toString();
            }
        }

        for (int likelihoodIndex = 0; likelihoodIndex < likelihoods.length; likelihoodIndex++) {
            if (genotypeLikelihoods[likelihoodIndex] > nextBestLikelihood && genotypeLikelihoods[likelihoodIndex] < bestLikelihood) {
                nextBestLikelihood = genotypeLikelihoods[likelihoodIndex];
            }
        }

        for (int likelihoodIndex = 0; likelihoodIndex < likelihoods.length; likelihoodIndex++) {
            if (refBase == Genotype.values()[likelihoodIndex].toString().charAt(0) &&
                refBase == Genotype.values()[likelihoodIndex].toString().charAt(1)) {
                refLikelihood = genotypeLikelihoods[likelihoodIndex];
            }
        }

        this.bestGenotype = bestGenotype;
        this.lodBtr = (float) (bestLikelihood - refLikelihood);
        this.lodBtnb = (float) (bestLikelihood - nextBestLikelihood);
    }

    public String getRefBasesFWD() {
        char[] b = { getReferenceBase() };
        return new String( b );
    }

    public char getRefSnpFWD() throws IllegalStateException { return getReferenceBase(); }
    public String getAltBasesFWD() { return getBestGenotype(); }
    public char getAltSnpFWD() throws IllegalStateException {
        String bases = getBestGenotype();
        if ( bases.charAt(0) != getRefSnpFWD() )
            return bases.charAt(0);
        else
            return bases.charAt(1);

    }
    public boolean isReference() { return ! isSNP(); }
    public boolean isSNP() { return getLodBtr() > 5; }
    public boolean isInsertion() { return false; }
    public boolean isDeletion() { return false; }
    public boolean isIndel() { return false; }
    public double getMAF() { return 0; }
    public double getHeterozygosity() { return 0; }
    public boolean isGenotype() { return true; }
    public double getVariationConfidence() { return getLodBtr(); }
    public double getConsensusConfidence() { return getLodBtnb(); }
    public List<String> getGenotype() throws IllegalStateException {
        return Arrays.asList(getBestGenotype());
    }
     
    public int getPloidy() throws IllegalStateException { return 2; }
    public boolean isBiallelic() { return true; }
}
