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
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers.CNV;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.PrintStream;
import java.util.List;

/**
 * Walks along reference and calculates the percent overlap with the BED file intervals for each -L interval.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = OverlapWithBedInIntervalWalker.INTERVALS_ROD_NAME, type = ReferenceOrderedDatum.class)})

public class OverlapWithBedInIntervalWalker extends RodWalker<CumulativeBaseOverlapCount, CumulativeBaseOverlapCount> {
    @Output
    protected PrintStream out;

    public final static String INTERVALS_ROD_NAME = "intervals";


    public boolean isReduceByInterval() {
        return true;
    }

    public void initialize() {
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public CumulativeBaseOverlapCount reduceInit() {
        return new CumulativeBaseOverlapCount();
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public CumulativeBaseOverlapCount map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        return new CumulativeBaseOverlapCount().addIntervals(tracker.getGATKFeatureMetaData(INTERVALS_ROD_NAME, true));
    }

    public CumulativeBaseOverlapCount reduce(CumulativeBaseOverlapCount add, CumulativeBaseOverlapCount runningCount) {
        if (add == null)
            add = new CumulativeBaseOverlapCount();

        return runningCount.addIn(add);
    }

    /**
     * @param results the genes found in each interval.
     */
    public void onTraversalDone(List<Pair<GenomeLoc, CumulativeBaseOverlapCount>> results) {
        for (Pair<GenomeLoc, CumulativeBaseOverlapCount> result : results ) {
            GenomeLoc loc = result.getFirst();

            CumulativeBaseOverlapCount overlapCount = result.getSecond();
            double meanOverlap = ((double) overlapCount.totalOverlapCount) / loc.size();

            out.println(loc + "\t" + meanOverlap);
        }
    }
}

class CumulativeBaseOverlapCount {
    public int totalOverlapCount;

    public CumulativeBaseOverlapCount() {
        this.totalOverlapCount = 0;
    }

    public CumulativeBaseOverlapCount addIn(CumulativeBaseOverlapCount other) {
        this.totalOverlapCount += other.totalOverlapCount;

        return this;
    }

    public CumulativeBaseOverlapCount addIntervals(List<GATKFeature> interval) {
        totalOverlapCount += interval.size();

        return this;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        sb.append(totalOverlapCount);

        return sb.toString();
    }
}