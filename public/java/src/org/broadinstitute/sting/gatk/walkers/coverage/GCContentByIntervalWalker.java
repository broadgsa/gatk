/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.PrintStream;
import java.util.List;

/**
 * Walks along reference and calculates the GC content for each interval.
 *
 *
 * <h2>Input</h2>
 * <p>
 *  A reference file
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 *  GC content calculations per interval.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T GCContentByInterval \
 *   -o output.txt \
 *   -L input.intervals
 * </pre>
 *
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})
@By(DataSource.REFERENCE)
public class GCContentByIntervalWalker extends LocusWalker<Long, Long> {
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

    public Long reduceInit() {
        return 0L;
    }

    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;
        int baseIndex = ref.getBaseIndex();
        return (baseIndex == BaseUtils.gIndex || baseIndex == BaseUtils.cIndex) ? 1L : 0L;
    }

    public Long reduce(Long toAdd, Long runningCount) {
        return runningCount + toAdd;
    }

    public void onTraversalDone(List<Pair<GenomeLoc, Long>> results) {
        for (Pair<GenomeLoc, Long> result : results ) {
            GenomeLoc loc = result.getFirst();
            Long gcCount = result.getSecond();

            double gcContent = (double) gcCount / loc.size();
            out.println(loc + "\t" + gcContent);
        }
    }
}
