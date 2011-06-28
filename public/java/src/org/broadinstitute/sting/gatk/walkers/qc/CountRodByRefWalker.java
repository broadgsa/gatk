/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.commandline.Argument;

/**
 * Prints out counts of the number of reference ordered data objects are
 * each locus for debugging RefWalkers.
 */
public class CountRodByRefWalker extends RefWalker<CountRodWalker.Datum, Pair<ExpandingArrayList<Long>, Long>> {
    @Argument(fullName = "verbose", shortName = "v", doc="If true, Countrod will print out detailed information about the rods it finds and locations", required = false)
    public boolean verbose = false;

    @Argument(fullName = "showSkipped", shortName = "s", doc="If true, CountRod will print out the skippped locations", required = false)
    public boolean showSkipped = false;

    CountRodWalker crw = new CountRodWalker();

    public void initialize() {
        crw.verbose = verbose;
        crw.showSkipped = showSkipped;
    }

    public CountRodWalker.Datum map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return crw.map(tracker, ref, context);
    }

    public Pair<ExpandingArrayList<Long>, Long> reduceInit() {
        return crw.reduceInit();
    }

    public Pair<ExpandingArrayList<Long>, Long> reduce(CountRodWalker.Datum point, Pair<ExpandingArrayList<Long>, Long> sum) {
        return crw.reduce(point, sum);
    }
}