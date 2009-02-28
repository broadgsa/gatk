/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.directed;

import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.picard.util.Interval;
import edu.mit.broad.picard.io.IoUtil;

import java.util.List;
import java.util.BitSet;
import java.io.File;

/**
 * Create a GenomeMask from an IntervalList or a file containing an IntervalList
 */
public class GenomeMaskFactory {

    public GenomeMask makeGenomeMaskFromIntervalList(IntervalList intervalList) {
        if (intervalList.getHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            intervalList.sort();
        }
        List<Interval> uniqueIntervals = intervalList.getUniqueIntervals();
        GenomeMask ret = new GenomeMask();

        SAMFileHeader samHeader = intervalList.getHeader();

        for (Interval interval : uniqueIntervals) {
            // TODO: Maybe figure out more intelligently how big the bitset might be?
            BitSet bitSet = ret.getOrCreate(samHeader.getSequenceIndex(interval.getSequence()), interval.getEnd() + 1);
            bitSet.set(interval.getStart(), interval.getEnd() + 1);
        }
        return ret;
    }

    public GenomeMask makeGenomeMaskFromIntervalList(File intervalListFile) {
        IoUtil.assertFileIsReadable(intervalListFile);
        IntervalList intervalList = IntervalList.fromFile(intervalListFile);
        return makeGenomeMaskFromIntervalList(intervalList);
    }
}
