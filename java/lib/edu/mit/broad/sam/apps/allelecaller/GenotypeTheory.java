package edu.mit.broad.sam.apps.allelecaller;

/**
 * Datastructure to hold a single genotype along with a likelihood.
 */
public class GenotypeTheory implements Comparable<GenotypeTheory> {
    private DiploidGenotype genotype;
    private double likelihood;

    public GenotypeTheory(final DiploidGenotype genotype, final double likelihood) {
        this.genotype = genotype;
        this.likelihood = likelihood;
    }

    public DiploidGenotype getGenotype() {
        return genotype;
    }

    public void setGenotype(final DiploidGenotype genotype) {
        this.genotype = genotype;
    }

    public double getLikelihood() {
        return likelihood;
    }

    public void setLikelihood(final double likelihood) {
        this.likelihood = likelihood;
    }

    /**
     * Genotype Theories are sorted first by descending likelihood (ie
     * the GenotypeTheory with biggest likelihood comes first).  Ties are
     * broken by lexical sorting of the genotypes themselves
     *
     */
    public int compareTo(final GenotypeTheory other) {
        if (this.getLikelihood() == other.getLikelihood()) {
            return this.getGenotype().compareTo(other.getGenotype());
        } else if (this.getLikelihood() > other.getLikelihood()) {
            return -1;
        } else {
            return 1;
        }
    }
}
