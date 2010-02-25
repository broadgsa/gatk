package org.broadinstitute.sting.gatk.refdata.utils;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class RODRecordList
 *         <p/>
 *         make the RODRecord list an interface, so we can stub in other implementations
 *         during testing.
 */
public interface RODRecordList extends List<ReferenceOrderedDatum>, Comparable<RODRecordList> {
    public GenomeLoc getLocation();
    public String getName();
}
