package org.broadinstitute.sting.utils.genotype;

/**
 * @author ebanks
 * Interface IDBacked
 *
 * this interface indicates that the genotype is
 * backed up by sample information.
 */
public interface IDBacked {

    /**
     *
     * @return returns the ID for this genotype
     */
    public String getID();

    /**
     *
     * @param   id    the ID for this genotype
     */
    public void setID(String id);

}