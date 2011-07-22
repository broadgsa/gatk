/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordCoordinateComparator;

/**
 * Extends Picard's Comparator for sorting SAMRecords by coordinate.  This one actually deals with unmapped reads
 * (among other things) sitting at the same position as their mates (so that they both can be put into the same set).
 */
public class SAMRecordCoordinateComparatorWithUnmappedReads extends SAMRecordCoordinateComparator {
    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int cmp = fileOrderCompare(samRecord1, samRecord2);
        if ( cmp != 0 )
            return cmp;

        // deal with unmapped reads
        if ( samRecord1.getReadUnmappedFlag() != samRecord2.getReadUnmappedFlag() )
            return (samRecord1.getReadUnmappedFlag()? 1: -1);        

        if ( samRecord1.getReadNegativeStrandFlag() != samRecord2.getReadNegativeStrandFlag() )
            return (samRecord1.getReadNegativeStrandFlag()? 1: -1);

        // even the names can be the same
        cmp = samRecord1.getReadName().compareTo(samRecord2.getReadName());
        if ( cmp != 0 )
            return cmp;

        if ( samRecord1.getDuplicateReadFlag() != samRecord2.getDuplicateReadFlag() )
            return (samRecord1.getDuplicateReadFlag()? -1: 1);

        if ( samRecord1.getReadPairedFlag() && samRecord2.getReadPairedFlag() && samRecord1.getFirstOfPairFlag() != samRecord2.getFirstOfPairFlag() )
            return (samRecord1.getFirstOfPairFlag()? -1: 1);

        // such a case was actually observed
        return samRecord1.getMappingQuality() - samRecord2.getMappingQuality();
    }
}
