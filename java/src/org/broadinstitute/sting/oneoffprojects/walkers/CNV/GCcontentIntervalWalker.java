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

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.*;

/**
 * Walks along all variant ROD loci and verifies the phasing from the reads for user-defined pairs of sites.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})

public class GCcontentIntervalWalker extends RodWalker<GCcounter, GCcounter> {
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

        List<GATKFeature> interval = tracker.getGATKFeatureMetaData("intervals", true);
        if (interval.size() != 1) {
            String error = "At " + ref.getLocus() + " : Must provide a track named 'intervals' with exactly ONE interval per locus in -L argument!";
            if (interval.size() < 1)
                throw new UserException(error);
            else // interval.size() > 1
                logger.warn(error);
        }
        GenomeLoc curInterval = interval.get(0).getLocation();

        GCcounter counter = new GCcounter();
        counter.calculateGCandAddIn(ref);
        counter.loc = curInterval;

        return counter;
    }

    public GCcounter reduce(GCcounter add, GCcounter runningCount) {
        if (add == null)
            add = new GCcounter();

        return runningCount.addIn(add);
    }

    /**
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(GCcounter result) {
        if (result.loc == null)
            return;

        double gcContent = (double) result.GCcount / result.totalCount;
        out.println(result.loc + "\t" + gcContent + "\t" + result.loc.size());
    }
}

class GCcounter {
    public int totalCount;
    public int GCcount;
    public GenomeLoc loc;

    public GCcounter() {
        this.totalCount = 0;
        this.GCcount = 0;
        this.loc = null;
    }

    public GCcounter addIn(GCcounter other) {
        this.totalCount += other.totalCount;
        this.GCcount += other.GCcount;

        if (other.loc != null && this.loc == null)
            this.loc = other.loc;

        return this;
    }

    public void calculateGCandAddIn(ReferenceContext ref) {
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
    }
}

