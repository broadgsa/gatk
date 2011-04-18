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

package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.PrintStream;

/**
 * Calculates basic statistics about the reference sequence itself
 */
public class FastaStatsWalker extends RefWalker<Byte, FastaStatsWalker.FastaStats> {
    @Output PrintStream out;

    protected class FastaStats {
        long nBases = 0, nRegBases = 0;
    }

    @Override
	public Byte map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        return ref.getBase();
	}

    @Override
    public FastaStats reduceInit() {
        return new FastaStats();
    }

    @Override
	public FastaStats reduce(Byte base, FastaStats stats) {
        stats.nBases++;
        if (BaseUtils.isRegularBase(base)) stats.nRegBases++;
        return stats;
	}

    @Override
    public void onTraversalDone(FastaStats sum) {
        out.printf("Total bases   %d%n", sum.nBases);
        out.printf("Regular bases %d%n", sum.nRegBases);
    }
}