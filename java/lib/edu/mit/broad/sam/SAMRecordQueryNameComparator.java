/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

/**
 * For "queryname" ordering of SAMRecords
 */
public class SAMRecordQueryNameComparator implements SAMRecordComparator {

    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        final int cmp = fileOrderCompare(samRecord1, samRecord2);
        if (cmp != 0) {
            return cmp;
        }
        if (samRecord1.getReadNegativeStrandFlag() == samRecord2.getReadNegativeStrandFlag()) {
            return 0;
        }
        return (samRecord1.getReadNegativeStrandFlag()? 1: -1);
    }

    /**
     * Less stringent compare method than the regular compare.  If the two records
     * are equal enough that their ordering in a sorted SAM file would be arbitrary,
     * this method returns 0.
     *
     * @return
     */
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        return samRecord1.getReadName().compareTo(samRecord2.getReadName());
    }
}
