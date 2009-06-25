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

package org.broadinstitute.sting.playground.gatk.walkers.indels;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
/**
 * Merges intervals based on reads which overlap them.
 */
@WalkerName("IntervalMerger")
@Requires({DataSource.READS})
public class IntervalMergerWalker extends ReadWalker<Integer,Integer> {

    @Argument(fullName="intervalsToMerge", shortName="intervals", doc="Intervals to merge", required=true)
    List<String> intervalsSource = null;
    @Argument(fullName="allow454Reads", shortName="454", doc="process 454 reads", required=false)
    public boolean allow454 = false;
    @Argument(fullName="maxIntervalSize", shortName="maxInterval", doc="max interval size", required=false)
    public int maxIntervalSize = 500;

    private LinkedList<GenomeLoc> intervals;
    private GenomeLoc firstInterval = null;

    @Override
    public void initialize() {
        intervals = parseIntervals(intervalsSource);
        firstInterval = intervals.removeFirst();
    }

    @Override
    public Integer map(char[] ref, SAMRecord read) {
        if ( firstInterval == null ||
             (!allow454 && Utils.is454Read(read, getToolkit().getEngine().getSAMHeader())) )
            return 0;

        GenomeLoc loc = GenomeLocParser.createGenomeLoc(read);

        // emit all intervals which we've passed
        while ( loc.isPast(firstInterval) ) {
            if ( firstInterval.getStop() - firstInterval.getStart() < maxIntervalSize)
                out.println(firstInterval);
            if ( intervals.size() > 0 ) {
                firstInterval = intervals.removeFirst();
            } else {
                firstInterval = null;
                return 0;
            }
        }

        // if we're not yet in the first interval, we're done
        if ( !loc.overlapsP(firstInterval))
            return 0;

        // at this point, we're in the first interval.
        // now we can merge any other intervals which we overlap
        while ( intervals.size() > 0 && loc.overlapsP(intervals.getFirst()) )
            firstInterval.setStop(intervals.removeFirst().getStop());

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
        if ( firstInterval != null &&
             firstInterval.getStop() - firstInterval.getStart() < maxIntervalSize)
            out.println(firstInterval);
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