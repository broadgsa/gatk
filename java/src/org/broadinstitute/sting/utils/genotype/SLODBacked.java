package org.broadinstitute.sting.utils.genotype;

/**
 * @author ebanks
 * Interface SLODBacked
 *
 * this interface indicates that the genotype is
 * backed up by strand lod information.
 */
public interface SLODBacked {

    /**
     *
     * @return returns the strand lod for this genotype or null if not set
     */
    public Double getSLOD();

    /**
     *
     * @param   slod    the strand lod for this genotype
     */
    public void setSLOD(double slod);

}