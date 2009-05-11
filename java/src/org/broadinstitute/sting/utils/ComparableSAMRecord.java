package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMRecord;

public class ComparableSAMRecord implements Comparable<ComparableSAMRecord> {

    private SAMRecord record;

    public ComparableSAMRecord(SAMRecord record) {
        this.record = record;
    }

    public SAMRecord getRecord() {
        return record;
    }

    public int compareTo(ComparableSAMRecord o) {
        // first sort by start position
        GenomeLoc myLoc = new GenomeLoc(record);
        GenomeLoc hisLoc = new GenomeLoc(o.getRecord());
        int comparison = myLoc.compareTo(hisLoc);
        // if the reads have the same start position, we must give a non-zero comparison
        // (because java Sets often require "consistency with equals")
        if ( comparison == 0 )
            comparison = record.getReadName().compareTo(o.getRecord().getReadName());
        return comparison;
    }
}