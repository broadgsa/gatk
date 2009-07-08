package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMRecord;

public class ComparableSAMRecord implements Comparable<ComparableSAMRecord> {

    private SAMRecord record;
    private GenomeLoc loc;

    public ComparableSAMRecord(SAMRecord record) {
        this.record = record;
        this.loc = new GenomeLoc(record);
    }

    public SAMRecord getRecord() {
        return record;
    }

    public int compareTo(ComparableSAMRecord o) {
        // first sort by start position
        int comparison = loc.compareTo(o.loc);
        // if the reads have the same start position, we must give a non-zero comparison
        // (because java Sets often require "consistency with equals")
        if ( comparison == 0 )
            comparison = record.getReadName().compareTo(o.getRecord().getReadName());
        // if the read names are the same, use the first of the pair if appropriate
        if ( comparison == 0 && record.getReadPairedFlag() )
            comparison = ( record.getFirstOfPairFlag() ? -1 : 1);
        return comparison;
    }

    public boolean equals(Object obj) {
        if ( !(obj instanceof ComparableSAMRecord) )
            return false;
        if ( this == obj )
            return true;

        ComparableSAMRecord csr = (ComparableSAMRecord)obj;
        if ( loc.compareTo(csr.loc) != 0 )
            return false;
        if ( !record.getReadName().equals(csr.getRecord().getReadName()) )
            return false;
        return ( record.getFirstOfPairFlag() == csr.record.getFirstOfPairFlag() );
    }
}