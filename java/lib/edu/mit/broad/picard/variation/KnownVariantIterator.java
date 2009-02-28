/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.variation;

import edu.mit.broad.sam.SAMSequenceRecord;

import java.util.Iterator;
import java.util.List;

/**
 * API for iterating over records representing known variations
 *
 * @author Kathleen Tibbetts
 */
public interface KnownVariantIterator extends Iterable<KnownVariant>, Iterator<KnownVariant>
{
    /**
     * Return the list of sequence dictionary (list of SAMSequenceRecords in order)
     * for this KnownVariantIterator
     *
     * @return The SAMSequenceRecords that comprise the sequence dictionary for this iterator, in order
     */
    public List<SAMSequenceRecord> getSequenceDictionary();    
}
