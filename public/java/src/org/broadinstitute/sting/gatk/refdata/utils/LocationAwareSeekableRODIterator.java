package org.broadinstitute.sting.gatk.refdata.utils;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * @author aaron
 *         <p/>
 *         Interface LocationAwareSeekableRODIterator
 *         <p/>
 *         combine iteration with a position aware interface
 */
public interface LocationAwareSeekableRODIterator extends CloseableIterator<RODRecordList> {
    public Object getHeader();

    public SAMSequenceDictionary getSequenceDictionary();

    public GenomeLoc peekNextLocation();

    public GenomeLoc position();

    public RODRecordList seekForward(GenomeLoc interval);    
}
