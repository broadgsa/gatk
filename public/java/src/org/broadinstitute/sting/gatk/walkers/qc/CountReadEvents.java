package org.broadinstitute.sting.gatk.walkers.qc;

import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Walks over the input data set, counting the number of read events (from the CIGAR operator)
 *
 * <h2>Input</h2>
 * <p>
 * One or more BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Number of reads events for each category
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountReadEvents \
 *   -o output.grp \
 *   -I input.bam \
 *   [-L input.intervals]
 * </pre>
 */

@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountReadEvents extends ReadWalker<Map<CigarOperator, ArrayList<Integer>> , Map<Integer, Map<CigarOperator, Long>>> {
    @Output (doc = "GATKReport table output")
    PrintStream out;

    public Map<CigarOperator, ArrayList<Integer>> map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
        return ReadUtils.getCigarOperatorForAllBases(read);
    }

    public Map<Integer, Map<CigarOperator, Long>> reduceInit() {
        return new HashMap<Integer, Map<CigarOperator, Long>>();
    }

    public Map<Integer, Map<CigarOperator, Long>> reduce(Map<CigarOperator, ArrayList<Integer>> value, Map<Integer, Map<CigarOperator, Long>> sum) {
        for (Map.Entry<CigarOperator, ArrayList<Integer>> entry : value.entrySet()) {
            CigarOperator op = entry.getKey();
            ArrayList<Integer> positions = entry.getValue();

            for (int p : positions) {
                Map<CigarOperator, Long> operatorCount = sum.get(p);
                if (operatorCount == null) {
                    operatorCount = new HashMap<CigarOperator, Long>();
                    sum.put(p, operatorCount);
                }

                Long count = operatorCount.get(op);
                if (count == null)
                    count = 0L;
                count++;
                operatorCount.put(op, count);
            }
        }
        return sum;
    }

    @Override
    public void onTraversalDone(Map<Integer, Map<CigarOperator, Long>> result) {
        GATKReport report = GATKReport.newSimpleReport("Events", "Position", "Event", "Observations");
        for (Map.Entry<Integer, Map<CigarOperator, Long>> entry : result.entrySet()) {
            int position = entry.getKey();
            Map<CigarOperator, Long> operatorCount = entry.getValue();

            for (Map.Entry<CigarOperator, Long> subEntry: operatorCount.entrySet()) {
                String operator = subEntry.getKey().name();
                Long observations = subEntry.getValue();
                report.addRow(position, operator, observations);
            }
        }
        report.print(out);
    }
}
