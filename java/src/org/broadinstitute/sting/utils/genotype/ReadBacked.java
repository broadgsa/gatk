package org.broadinstitute.sting.utils.genotype;

import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * @author aaron
 *
 * Interface ReadBacked
 *
 * This interface indicates that a genotype or variant is backed by reads
 */
public interface ReadBacked {
    /**
     * get the SAM records that back this genotype call
     * @return a list of SAM records
     */
    public List<SAMRecord> getReads();

    /**
     * get the count of reads
     * @return the number of reads we're backed by
     */
    public int getReadCount();
}
