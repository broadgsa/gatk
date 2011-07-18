/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.sam;

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
        // first sort by start position -- with not coverflow because both are guaranteed to be positive.
        int comparison = record.getAlignmentStart() - o.record.getAlignmentStart();
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
        if(record.getAlignmentStart() != csr.record.getAlignmentStart())
            return false;
        if ( !record.getReadName().equals(csr.getRecord().getReadName()) )
            return false;
        return ( record.getFirstOfPairFlag() == csr.record.getFirstOfPairFlag() );
    }
}