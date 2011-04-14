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
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.*;

/**
 * Walks along reference and calculates the GC content for each interval.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})

@By(DataSource.REFERENCE)

public class GCcontentIntervalWalker extends LocusWalker<GCcounter, GCcounter> {
    @Output
    protected PrintStream out;

    public boolean isReduceByInterval() {
        return true;
    }

    public void initialize() {
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public GCcounter reduceInit() {
        return new GCcounter();
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public GCcounter map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        return new GCcounter().calculateGCandAddIn(ref);
    }

    public GCcounter reduce(GCcounter add, GCcounter runningCount) {
        if (add == null)
            add = new GCcounter();

        return runningCount.addIn(add);
    }

    /**
     * @param results the GC content observed for each interval.
     */
    public void onTraversalDone(List<Pair<GenomeLoc, GCcounter>> results) {
        for (Pair<GenomeLoc, GCcounter> result : results ) {
            GenomeLoc loc = result.getFirst();
            GCcounter counter = result.getSecond();

            double gcContent = (double) counter.GCcount / counter.totalCount;
            out.println(loc + "\t" + gcContent + "\t" + loc.size());
        }
    }
}

class GCcounter {
    public int totalCount;
    public int GCcount;

    public GCcounter() {
        this.totalCount = 0;
        this.GCcount = 0;
    }

    public GCcounter addIn(GCcounter other) {
        this.totalCount += other.totalCount;
        this.GCcount += other.GCcount;

        return this;
    }

    public GCcounter calculateGCandAddIn(ReferenceContext ref) {
        for (byte base : ref.getBases()) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);

            boolean baseIsGC = (baseIndex == BaseUtils.gIndex || baseIndex == BaseUtils.cIndex);
            boolean baseIsAT = (baseIndex == BaseUtils.aIndex || baseIndex == BaseUtils.tIndex);
            if (baseIsGC || baseIsAT) {
                totalCount++;
                if (baseIsGC)
                    GCcount++;
            }
        }

        return this;
    }
}