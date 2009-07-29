package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.Pair;

/**
 * @author aaron
 *         <p/>
 *         Interface GenotypeCall
 *         <p/>
 *         Genotype call interface, for indicating that a genotype is
 *         also a genotype call.
 */
public interface GenotypeCall extends GenotypeLocus{
    /**
     * get the confidence
     *
     * @return a ConfidenceScore representing the LOD score that this genotype was called with
     */
    public ConfidenceScore getConfidence();

    /**
     * gets the reference base
     *
     * @return the reference base we represent
     */
    public char getReferencebase();

    /**
     * get the best vrs the next best genotype LOD score
     * @return the genotype, and a LOD for best - next
     */
    public Pair<Genotype,ConfidenceScore> getBestVrsNext();

    /**
     * get the best vrs the reference allele.
     * @return the genotype, and a LOD for best - ref.  The best may be ref, unless you've checked
     * with is variation
     */
    public Pair<Genotype,ConfidenceScore> getBestVrsRef();

    /**
     * check to see if this call is a variant, i.e. not homozygous reference
     * @return true if it's not hom ref, false otherwise
     */
    public boolean isVariation();

    /**
     * return genotype locus, with our data
     */
    public GenotypeLocus toGenotypeLocus();
}

