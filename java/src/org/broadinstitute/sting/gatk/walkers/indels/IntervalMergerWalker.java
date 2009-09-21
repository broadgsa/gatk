/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.filters.Platform454Filter;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;


/**
 * Merges intervals based on reads which overlap them.
 */
@WalkerName("IntervalMerger")
@Requires({DataSource.READS, DataSource.REFERENCE})
@ReadFilters({Platform454Filter.class, ZeroMappingQualityReadFilter.class})
public class IntervalMergerWalker extends ReadWalker<Integer,Integer> {

    @Argument(fullName="intervalsToMerge", shortName="intervals", doc="Intervals to merge", required=true)
    List<String> intervalsSource = null;
    @Argument(fullName="allow454Reads", shortName="454", doc="process 454 reads", required=false)
    boolean allow454 = false;
    @Argument(fullName="maxIntervalSize", shortName="maxInterval", doc="max interval size", required=false)
    int maxIntervalSize = 500;

    private LinkedList<GenomeLoc> intervals;
    private GenomeLoc currentInterval = null;
    private boolean currentIntervalIsUsed = false;

    @Override
    public void initialize() {
        intervals = parseIntervals(intervalsSource);
        currentInterval = (intervals.size() > 0 ? intervals.removeFirst() : null);
    }

    @Override
    public Integer map(char[] ref, SAMRecord read) {
        if ( currentInterval == null ||
             (!allow454 && Utils.is454Read(read)) )
            return 0;

        GenomeLoc loc = GenomeLocParser.createGenomeLoc(read);

        // ignore all intervals which we've passed
        while ( loc.isPast(currentInterval) ) {
            if ( currentIntervalIsUsed && currentInterval.getStop() - currentInterval.getStart() < maxIntervalSize) {
                out.println(currentInterval);
                currentIntervalIsUsed = false;
            }
            if ( intervals.size() > 0 ) {
                currentInterval = intervals.removeFirst();
            } else {
                currentInterval = null;
                return 0;
            }
        }

        // if we're not yet in the current interval, we're done
        if ( !loc.overlapsP(currentInterval))
            return 0;

        // at this point, we're in the current interval.
        // now we can merge any other intervals which we overlap
        currentIntervalIsUsed = true;
        while ( intervals.size() > 0 && loc.overlapsP(intervals.getFirst()) )
            currentInterval = GenomeLocParser.setStop(currentInterval, intervals.removeFirst().getStop());

        return 1;
    }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce( Integer value, Integer accum ) {
        return accum + value;
    }

    @Override
    public void onTraversalDone( Integer value ) {
        if ( currentInterval != null && currentIntervalIsUsed &&
             currentInterval.getStop() - currentInterval.getStart() < maxIntervalSize)
            out.println(currentInterval);
    }

    /**
     * Load the intervals directly from the command-line or from file, as appropriate.
     * Merge overlapping intervals.
     * @param intervalsSource Source of intervals.
     * @return a linked list of sorted, merged intervals.
     */
    private LinkedList<GenomeLoc> parseIntervals(List<String> intervalsSource) {
        List<GenomeLoc> parsedIntervals = GenomeAnalysisEngine.parseIntervalRegion(intervalsSource);
        GenomeLocSortedSet intervalSortedSet = new GenomeLocSortedSet();
        for ( GenomeLoc parsedInterval : parsedIntervals )
            intervalSortedSet.addRegion(parsedInterval);

        return new LinkedList<GenomeLoc>( intervalSortedSet );
    }
}
