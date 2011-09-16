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

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.PrintStream;

/**
 * Renders a new reference in FASTA format consisting of only those loci provided in the input data set.
 *
 * <p>
 * The output format can be partially controlled using the provided command-line arguments.
 * Specify intervals with the usual -L argument to output only the reference bases within your intervals.
 * Overlapping intervals are automatically merged; reference bases for each disjoint interval will be output as a
 * separate fasta sequence (named numerically in order).
 *
 * <h2>Input</h2>
 * <p>
 * The reference and requested intervals.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A fasta file representing the requested intervals.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T FastaReference \
 *   -o output.fasta \
 *   -L input.intervals
 * </pre>
 *
 */
@WalkerName("FastaReferenceMaker")
public class FastaReferenceWalker extends RefWalker<Pair<GenomeLoc, String>, GenomeLoc> {

    @Output PrintStream out;

    @Argument(fullName="lineWidth", shortName="lw", doc="Maximum length of sequence to write per line", required=false)
    public int fastaLineWidth=60;

    /**
     *  Please note that when using this argument adjacent intervals will automatically be merged.
     */
    @Argument(fullName="rawOnelineSeq", shortName="raw", doc="Print sequences with no FASTA header lines, one line per interval (i.e. lineWidth = infinity)", required=false)
    public boolean fastaRawSeqs=false;

    protected FastaSequence fasta;

    public void initialize() {
        if (fastaRawSeqs) fastaLineWidth = Integer.MAX_VALUE;
        fasta = new FastaSequence(out, fastaLineWidth, fastaRawSeqs);
    }

	public Pair<GenomeLoc, String> map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        return new Pair<GenomeLoc, String>(context.getLocation(), String.valueOf((char)ref.getBase()));
	}

    public GenomeLoc reduceInit() {
        return null;
    }

	public GenomeLoc reduce(Pair<GenomeLoc, String> value, GenomeLoc sum) {
        if ( value == null )
            return sum;

        // if there is no interval to the left, then this is the first one
        if ( sum == null ) {
            sum = value.first;
            fasta.append(value.second);
        }
        // if the intervals don't overlap, print out the leftmost one and start a new one
        // (end of contig or new interval)
        else if ( value.first.getStart() != sum.getStop() + 1 ) {
            fasta.flush();
            sum = value.first;
            fasta.append(value.second);
        }
        // otherwise, merge them
        else {
            sum = getToolkit().getGenomeLocParser().setStop(sum, value.first.getStop());
            fasta.append(value.second);
        }
		return sum;
	}

    public void onTraversalDone(GenomeLoc sum) {
        fasta.flush();
    }
}