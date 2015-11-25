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

package org.broadinstitute.gatk.tools.walkers.fasta;

import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RefWalker;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.PrintStream;

/**
 * Calculate basic statistics about the reference sequence itself
 *
 * <p>These are very basic statistics: total number of bases and number of "regular" bases (i.e. A, C, T or G).</p>
 *
 * <h3>Input</h3>
 * <p>
 * A FASTA reference file.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Base counts are written to file if an output file name is given (with -o), otherwise output to stdout.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T FastaStats \
 *   -R reference.fasta \
 *   [-o output.txt]
 * </pre>
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class FastaStats extends RefWalker<Byte, FastaStats.FastaStatistics> {
    @Output PrintStream out;

    protected class FastaStatistics {
        long nBases = 0, nRegBases = 0;
    }

    @Override
	public Byte map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        return ref.getBase();
	}

    @Override
    public FastaStatistics reduceInit() {
        return new FastaStatistics();
    }

    @Override
	public FastaStatistics reduce(Byte base, FastaStatistics stats) {
        stats.nBases++;
        if (BaseUtils.isRegularBase(base)) stats.nRegBases++;
        return stats;
	}

    @Override
    public void onTraversalDone(FastaStatistics sum) {
        out.printf("Total bases   %d%n", sum.nBases);
        out.printf("Regular bases %d%n", sum.nRegBases);
    }
}