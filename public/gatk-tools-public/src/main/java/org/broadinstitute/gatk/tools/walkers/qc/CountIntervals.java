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

package org.broadinstitute.gatk.tools.walkers.qc;

import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RefWalker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.PrintStream;
import java.util.Collections;
import java.util.List;

/**
 * Count contiguous regions in an interval list
 *
 * <p>When the GATK reads in intervals from an intervals list, any intervals that overlap each other get merged into
 * a single interval spanning the original ones. For example, if you have the following intervals:
 * <ul><li>
 *     20:1-2000
 * </li><li>
 *     20:1500-3000
 * </li></ul>
 * They will be merged into a single interval:
 * <ul><li>20:1-3000</li></ul>
 *
 * This tool allows you to check, for a given list of intervals, how many separate intervals the GATK will actually
 * distinguish at runtime.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more ROD files containing intervals to check.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Number of separate intervals identified by GATK after merging overlapping intervals.
 * </p>
 *
 * You can use the -numOverlaps argument to find out how many cases you have of a specific number of overlaps.
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CountIntervals \
 *   -R reference.fasta \
 *   -o output.txt \
 *   -check intervals.list
 * </pre>
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class CountIntervals extends RefWalker<Long, Long> {
    @Output
    PrintStream out;

    @Input(fullName="check", shortName = "check", doc="Any number of RODs", required=false)
    public List<RodBinding<Feature>> features = Collections.emptyList();

    @Argument(fullName="numOverlaps",shortName="no",doc="Count all occurrences of X or more overlapping intervals; defaults to 2", required=false)
    int numOverlaps = 2;

    public Long reduceInit() {
        return 0l;
    }

    public boolean isReduceByInterval() { return true; }

    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return null;
        }

        List<Feature> checkIntervals = tracker.getValues(features);
        return (long) checkIntervals.size();
    }

    public Long reduce(Long loc, Long prev) {
        if ( loc == null ) {
            return 0l;
        } else {
            return Math.max(prev,loc);
        }
    }

    public void onTraversalDone(List<Pair<GenomeLoc,Long>> finalReduce) {
        long count = 0;
        for ( Pair<GenomeLoc,Long> g : finalReduce ) {
            if ( g.second >= numOverlaps) {
                count ++;
            }
        }
        out.printf("Number of contiguous intervals: %d",count);
    }
}
