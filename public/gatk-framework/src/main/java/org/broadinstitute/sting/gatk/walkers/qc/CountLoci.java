/*
* Copyright (c) 2012 The Broad Institute
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

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.NanoSchedulable;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;

import java.io.PrintStream;

/**
 * Walks over the input data set, calculating the total number of covered loci for diagnostic purposes.
 *
 * <p>
 * This is the simplest example of a locus walker.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more BAM files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Number of loci traversed. If an output file name is provided, then the result will be written to that file.
 * Otherwise it will be sent to standard console output.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -T CountLoci \
 *   -R ref.fasta \
 *   -I input.bam \
 *   -o output.txt \
 *   [-L input.intervals]
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class CountLoci extends LocusWalker<Integer, Long> implements TreeReducible<Long>, NanoSchedulable {
    @Output
    PrintStream out;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return 1;
    }

    public Long reduceInit() { return 0l; }

    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }

    /**
     * Reduces two subtrees together.  In this case, the implementation of the tree reduce
     * is exactly the same as the implementation of the single reduce.
     */
    public Long treeReduce(Long lhs, Long rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone( Long c ) {
        out.println(c);
    }
}
