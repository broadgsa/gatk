package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.NanoSchedulable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Walks over the input data set, calculating the number of reads seen for diagnostic purposes.
 *
 * <p>
 * Can also count the number of reads matching a given criterion using read filters (see the
 * --read-filter command line argument).  Simplest example of a read-backed analysis.
 *
 *
 * <h2>Input</h2>
 * <p>
 * One or more BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Number of reads seen.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountReads \
 *   -I input.bam \
 *   [-L input.intervals]
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountReads extends ReadWalker<Integer, Integer> implements NanoSchedulable {
    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
        return 1;
    }

    @Override public Integer reduceInit() { return 0; }
    @Override public Integer reduce(Integer value, Integer sum) { return value + sum; }
}
