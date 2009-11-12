package org.broadinstitute.sting.utils.genotype;

import java.util.Map;

/**
 * @author ebanks
 * Interface AlleleBalanceBacked
 *
 * this interface indicates that the genotype is
 * backed up by allele balance information.
 */
public interface ArbitraryFieldsBacked {

    /**
     *
     * @return returns the list of arbitrary fields to emit
     */
    public Map<String, String> getFields();

    /**
     *
     * @param   fields  a list of arbitrary fields to emit
     */
    public void setFields(Map<String, String> fields);
}