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
 * Comparator for sorting SAMRecords by coordinate.  Note that the header is required because
 * the order of sequences in the header defines the major sort order.
 */
public class SAMRecordCoordinateComparator implements SAMRecordComparator {
    private final SAMFileHeader header;
    public SAMRecordCoordinateComparator(final SAMFileHeader header) {
        this.header = header;
    }
    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        final int cmp = fileOrderCompare(samRecord1, samRecord2);
        if (cmp != 0) {
            return cmp;
        }
        if (samRecord1.getReadNegativeStrandFlag() == samRecord2.getReadNegativeStrandFlag()) {
            return samRecord1.getReadName().compareTo(samRecord2.getReadName());
        }
        else {
            return (samRecord1.getReadNegativeStrandFlag()? 1: -1);
        }



    }

    /**
     * Less stringent compare method than the regular compare.  If the two records
     * are equal enough that their ordering in a sorted SAM file would be arbitrary,
     * this method returns 0.  If read is paired and unmapped, use the mate mapping to sort.
     *
     * @return
     */
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int refIndex1 = samRecord1.getReferenceIndex(header);
        int refIndex2 = samRecord2.getReferenceIndex(header);
        if (refIndex1 == -1) {
            return (refIndex2 == -1? 0: 1);
        } else if (refIndex2 == -1) {
            return -1;
        }
        int cmp = refIndex1 - refIndex2;
        if (cmp != 0) {
            return cmp;
        }
        return samRecord1.getAlignmentStart() - samRecord2.getAlignmentStart();
    }
}
