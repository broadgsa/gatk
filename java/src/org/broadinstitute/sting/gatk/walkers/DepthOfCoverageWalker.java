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

package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;

/**
 * Display the depth of coverage at a given locus.
 */
public class DepthOfCoverageWalker extends LocusWalker<Integer, Pair<Long, Long>> {
    @Argument(fullName="suppressLocusPrinting",doc="Suppress printing",required=false)
    public boolean suppressPrinting = false;

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        if ( !suppressPrinting )
            out.printf("%s: %d%n", context.getLocation(), context.getReads().size() );
        return context.getReads().size();
    }

    public boolean isReduceByInterval() {
        return true;
    }

    public Pair<Long, Long> reduceInit() { return new Pair<Long,Long>(0l,0l); }

    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum) {
        long left = value.longValue() + sum.getFirst();
        long right = sum.getSecond() + 1l;
        return new Pair<Long,Long>(left, right);
    }

    public void onTraversalDone(Pair<Long, Long> result) {
        out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n", 
                ((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
    }
}