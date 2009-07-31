package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.confidence.ConfidenceScore;
import org.broadinstitute.sting.utils.genotype.variant.Variant;

/**
 * @author aaron
 *         <p/>
 *         Class GenotypeLikelihood
 *         <p/>
 *         This class emcompasses all the basic information about a genotype
 */
public interface Genotype {
    /**
     * get the confidence score
     *
     * @return get the confidence score that we're based on
     */
    public ConfidenceScore getConfidenceScore();

    /**
     * get the bases that represent this
     *
     * @return the bases, as a string
     */
    public String getBases();

    /**
     * get the ploidy
     *
     * @return the ploidy value
     */
    public int getPloidy();

    /**
     * Returns true if both observed alleles are the same (regardless of whether they are ref or alt)
     *
     * @return true if we're homozygous, false otherwise
     */
    public boolean isHom();

    /**
     * Returns true if observed alleles differ (regardless of whether they are ref or alt)
     *
     * @return true if we're het, false otherwise
     */
    public boolean isHet();

    /**
     * get the genotype's location
     *
     * @return a GenomeLoc representing the location
     */
    public GenomeLoc getLocation();

    /**
     * returns true if the genotype is a point genotype, false if it's a indel / deletion
     *
     * @return true is a SNP
     */
    public boolean isPointGenotype();

    /**
     * given the reference, are we a variant? (non-ref)
     *
     * @param ref the reference base or bases
     *
     * @return true if we're a variant
     */
    public boolean isVariant(char ref);

    /**
     * return this genotype as a variant
     *
     * @return
     */
    public Variant toVariant();
}
