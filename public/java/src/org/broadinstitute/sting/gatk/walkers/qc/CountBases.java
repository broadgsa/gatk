package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Walks over the input data set, calculating the number of bases seen for diagnostic purposes.
 *
 * <h2>Input</h2>
 * <p>
 * One or more BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Number of bases seen.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountBases \
 *   -o output.txt \
 *   -I input.bam \
 *   [-L input.intervals]
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountBases extends ReadWalker<Integer, Long> {
    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {

        return read.getReadLength();
    }

    public Long reduceInit() { return 0L; }

    public Long reduce(Integer value, Long sum) {
        return (long) value + sum;
    }
}
