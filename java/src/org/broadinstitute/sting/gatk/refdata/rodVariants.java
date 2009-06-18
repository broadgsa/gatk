package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.IOException;

public class rodVariants extends BasicReferenceOrderedDatum {
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
            loc = new GenomeLoc(parts[0], Long.valueOf(parts[1]));
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
        return String.format("%s %c %d %d %s %4.4f %4.4f %f %f %f %f %f %f %f %f %f %f",
                loc,
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
}
