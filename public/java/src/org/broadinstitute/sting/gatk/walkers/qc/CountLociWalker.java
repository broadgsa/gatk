package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;

import java.io.PrintStream;

/**
 * Walks over the input data set, calculating the total number of covered loci for diagnostic purposes.
 *
 * <p>
 * Simplest example of a locus walker.
 *
 *
 * <h2>Input</h2>
 * <p>
 * One or more BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Number of loci traversed.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountLoci \
 *   -o output.txt \
 *   -I input.bam \
 *   [-L input.intervals]
 * </pre>
 *
 */
public class CountLociWalker extends LocusWalker<Integer, Long> implements TreeReducible<Long> {
    @Output(doc="Write count to this file instead of STDOUT")
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
