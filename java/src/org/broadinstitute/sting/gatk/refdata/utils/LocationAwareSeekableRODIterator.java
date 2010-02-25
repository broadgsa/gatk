package org.broadinstitute.sting.gatk.refdata.utils;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Iterator;
import java.util.List;

/**
 * @author aaron
 *         <p/>
 *         Interface LocationAwareSeekableRODIterator
 *         <p/>
 *         combine iteration with a position aware interface
 */
public interface LocationAwareSeekableRODIterator extends Iterator<RODRecordList> {
    public GenomeLoc peekNextLocation();

    public GenomeLoc position();

    public RODRecordList seekForward(GenomeLoc interval);

}
