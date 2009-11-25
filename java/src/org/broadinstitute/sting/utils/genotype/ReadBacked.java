package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * @author ebanks
 *
 * Interface ReadBacked
 *
 * This interface indicates that a genotype or variant is backed by reads
 */
public interface ReadBacked {
    /**
     * get the pileup that backs this genotype call
     * @return a pileup
     */
    public ReadBackedPileup getPileup();

    /**
     *
     * @param   pileup    the pileup for this genotype
     */
    public void setPileup(ReadBackedPileup pileup);

    /**
     * get the count of reads
     * @return the number of reads we're backed by
     */
    public int getReadCount();
}
