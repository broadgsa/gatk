/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.coverage;

import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.PrintStream;
import java.util.List;

/**
 * Calculates the GC content of the reference sequence for each interval
 *
 *
 * <h3>Input</h3>
 * <p>
 *  A reference file
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 *  GC content calculations per interval.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T GCContentByInterval \
 *   -R reference.fasta \
 *   -o output.txt \
 *   -L input.intervals
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})
@By(DataSource.REFERENCE)
public class GCContentByInterval extends LocusWalker<Long, Long> {
    @Output
    protected PrintStream out;

    public boolean isReduceByInterval() {
        return true;
    }

    public void initialize() {
    }

    public Long reduceInit() {
        return 0L;
    }

    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;
        int baseIndex = ref.getBaseIndex();
        return (baseIndex == BaseUtils.Base.G.ordinal() || baseIndex == BaseUtils.Base.C.ordinal()) ? 1L : 0L;
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
