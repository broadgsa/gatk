/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.interval.IntervalFileMergingIterator;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.commandline.Argument;

import java.io.File;
import java.util.*;

@By(DataSource.READS)
// walker to count realigned reads
public class RealignedReadCounter extends ReadWalker<Integer, Integer> {

    public static final String ORIGINAL_CIGAR_TAG = "OC";
    public static final String ORIGINAL_POSITION_TAG = "OP";

    @Argument(fullName="targetIntervals", shortName="targetIntervals", doc="intervals file output from RealignerTargetCreator", required=true)
    protected String intervalsFile = null;

    // the intervals input by the user
    private Iterator<GenomeLoc> intervals = null;

    // the current interval in the list
    private GenomeLoc currentInterval = null;

    private long updatedIntervals = 0, updatedReads = 0;
    private boolean intervalWasUpdated = false;

    public void initialize() {
        // prepare to read intervals one-by-one, as needed (assuming they are sorted).
        intervals = new IntervalFileMergingIterator( new File(intervalsFile), IntervalMergingRule.OVERLAPPING_ONLY );
        currentInterval = intervals.hasNext() ? intervals.next() : null;
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        if ( currentInterval == null ) {
            return 0;
        }

        GenomeLoc readLoc = GenomeLocParser.createGenomeLoc(read);
        // hack to get around unmapped reads having screwy locations
        if ( readLoc.getStop() == 0 )
            readLoc = GenomeLocParser.createGenomeLoc(readLoc.getContigIndex(), readLoc.getStart(), readLoc.getStart());

        if ( readLoc.isBefore(currentInterval) || ReadUtils.is454Read(read) )
            return 0;

        if ( readLoc.overlapsP(currentInterval) ) {
            if ( doNotTryToClean(read) )
                return 0;

            if ( read.getAttribute(ORIGINAL_CIGAR_TAG) != null ) {
                String newCigar = (String)read.getAttribute(ORIGINAL_CIGAR_TAG);
                // deal with an old bug
                if ( read.getCigar().toString().equals(newCigar) ) {
                    //System.out.println(currentInterval + ": " + read.getReadName() + " " + read.getCigarString() + " " + newCigar);
                    return 0;
                }

                if ( !intervalWasUpdated ) {
                    intervalWasUpdated = true;
                    updatedIntervals++;
                }
                updatedReads++;

            }
        } else {
            do {
                intervalWasUpdated = false;
                currentInterval = intervals.hasNext() ? intervals.next() : null;
            } while ( currentInterval != null && currentInterval.isBefore(readLoc) );
        }

        return 0;
    }

    private boolean doNotTryToClean(SAMRecord read) {
        return read.getReadUnmappedFlag() ||
                read.getNotPrimaryAlignmentFlag() ||
                read.getReadFailsVendorQualityCheckFlag() ||
                read.getMappingQuality() == 0 ||
                read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ||
                (BadMateFilter.hasBadMate(read));
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        System.out.println(updatedIntervals + " intervals were updated");
        System.out.println(updatedReads + " reads were updated");
    }
}
