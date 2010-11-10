/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broad.tribble.util.variantcontext.VariantContext;

import java.util.List;
import java.util.Collection;
import java.util.ArrayList;
import java.io.PrintStream;
import java.io.InputStream;
import java.io.File;
import java.io.FileNotFoundException;

import net.sf.picard.util.IntervalList;
import net.sf.picard.util.Interval;

/**
 * Counts the number of contiguous regions the walker traverses over. Slower than it needs to be, but
 * very useful since overlapping intervals get merged, so you can count the number of intervals the GATK merges down to.
 * This was its very first use.
 */
public class MarkIntervals extends RodWalker<Long, Long> {
    @Output
    File out;

    @Input(doc="List of bad SNP sites", required=true)
    File locs;

    List<GenomeLoc> badSites = new ArrayList<GenomeLoc>();
    IntervalList intervalList;

    @Argument(doc="Should we match intervals or just the sites?", required=false)
    boolean matchJustSites = false;

    public void initialize() {
        if ( this.getToolkit().getArguments().intervals.size() != 1 )
            throw new UserException("This walker only works with a single -L argument provided, sorry");

        File intervalFile = new File(this.getToolkit().getArguments().intervals.get(0));
        try {
            intervalList = IntervalList.fromFile(intervalFile);
        } catch ( Exception ex ) {
            throw new UserException.MalformedFile(intervalFile, "Couldn't read interval file", ex);
        }

        try {
            for ( String line : new XReadLines(locs, true) ) {
                String parts[] = line.split(":");
                badSites.add(getToolkit().getGenomeLocParser().createGenomeLoc(parts[0], Integer.valueOf(parts[1])));
            }
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(locs, e);
        }
    }

    public Long reduceInit() {
        return 0l;
    }

    public boolean isReduceByInterval() { return true; }

    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return null;
        }

        return 0L;
    }

    public Long reduce(Long loc, Long prev) {
        return 0L;
    }

    public void onTraversalDone(List<Pair<GenomeLoc,Long>> finalReduce) {
        long nBadIntervals = 0;
        long nBadBases = 0;

        IntervalList badIntervals = new IntervalList(intervalList.getHeader());

        if ( matchJustSites ) {
            for ( GenomeLoc loc : badSites ) {
                nBadIntervals++;
                nBadBases += loc.size();
                badIntervals.add(new Interval(loc.getContig(), (int)loc.getStart(), (int)loc.getStop(), false, "nBadSites_" + 1));
            }
        } else {
            for ( Pair<GenomeLoc,Long> g : finalReduce ) {
                int overlaps = 0;
                GenomeLoc interval = g.getFirst();
                for ( GenomeLoc loc : badSites ) {
                    if ( interval.overlapsP(loc) ) {
                        logger.info(String.format("Overlapping %s with bad site %s", interval, loc));
                        overlaps++;
                    }
                }

                if ( overlaps > 0 ) {
                    nBadIntervals++;
                    nBadBases += interval.size();
                    badIntervals.add(new Interval(interval.getContig(), (int)interval.getStart(), (int)interval.getStop(), false, "nBadSites_" + overlaps));
                    //out.printf("%s %d%n", interval, overlaps);
                }
            }
        }

        logger.info("No. intervals marked as bad: " + nBadIntervals);
        logger.info("No.     bases marked as bad: " + nBadBases);
        badIntervals.write(out);
    }
}