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

package org.broadinstitute.sting.playground.gatk.walkers.diagnostics;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import net.sf.samtools.SAMRecord;

import java.util.*;
import java.io.PrintStream;

/**
 * Compute quality score distribution
 */
public class QualityScoreDistribution extends LocusWalker<Integer, Integer> {
    @Output
    PrintStream out;

    private HashMap<String, long[]> qualDists;

    public void initialize() {
        qualDists = new HashMap<String, long[]>();

        qualDists.put("all", new long[QualityUtils.MAX_QUAL_SCORE]);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        for (int i = 0; i < reads.size(); i++) {
            byte qual = reads.get(i).getBaseQualities()[offsets.get(i)];
            String name = reads.get(i).getReadGroup().getReadGroupId();

            if (!qualDists.containsKey(name)) {
                qualDists.put(name, new long[QualityUtils.MAX_QUAL_SCORE]);
            }

            qualDists.get(name)[qual]++;
            qualDists.get("all")[qual]++;
        }

        return null;
    }

    public Integer reduceInit() {
        return null;
    }

    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer result) {
        Set<String> names = qualDists.keySet();
        HashMap<String, Long> norms = new HashMap<String, Long>();

        for (String name : names) {
            long norm = 0;
            for (int qual = 0; qual < QualityUtils.MAX_QUAL_SCORE; qual++) {
                norm += qualDists.get(name)[qual];
            }

            norms.put(name, norm);
        }

        out.printf("Q");
        for (String name : names) {
            out.printf("\t%s", name);
        }
        out.println();

        for (int qual = 0; qual < QualityUtils.MAX_QUAL_SCORE; qual++) {
            out.printf("%d", qual);

            for (String name : names) {
                out.printf("\t%f", ((float) qualDists.get(name)[qual])/((float) norms.get(name)));
            }

            out.println();
        }
    }
}