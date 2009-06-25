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
 * User: ebanks
 * Date: Jun 10, 2009
 * Time: 2:40:19 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Merges intervals based on reads which overlap them.
 */
@WalkerName("IntervalMerger")
@Requires({DataSource.READS})
public class IntervalMergerWalker extends ReadWalker<Integer,Integer> {

    @Argument(fullName="intervalsToMerge", shortName="intervals", doc="Intervals to merge", required=true)
    String intervalsSource = null;
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
    private LinkedList<GenomeLoc> parseIntervals(String intervalsSource) {
        List<GenomeLoc> parsedIntervals = GenomeAnalysisEngine.parseIntervalRegion(intervalsSource,false);
        GenomeLocSortedSet intervalSortedSet = new GenomeLocSortedSet();
        for ( GenomeLoc parsedInterval : parsedIntervals )
            intervalSortedSet.addRegion(parsedInterval);

        return new LinkedList<GenomeLoc>( intervalSortedSet );
    }
}