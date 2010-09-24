package org.broadinstitute.sting.playground.gatk.walkers.diagnostics;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.PrintStream;
import java.util.List;

/**
 * Computes the read error rate per position in read (in the original 5'->3' orientation that the read had coming off the machine)
 */
public class ErrorRatePerCycle extends LocusWalker<Integer, Integer> {
    @Output PrintStream out;
    @Argument(fullName="min_base_quality_score", shortName="mbq", doc="Minimum base quality required to consider a base for calling (default: 0)", required=false) public Integer MIN_BASE_QUAL = 0;
    @Argument(fullName="min_mapping_quality_score", shortName="mmq", doc="Minimum read mapping quality required to consider a read for calling (default: 0)", required=false) public Integer MIN_MAPPING_QUAL = 0;

    private GATKReport report;
    private String reportName = "ErrorRatePerCycle";
    private String reportDescription = "The error rate per sequenced position in the reads";

    public void initialize() {
        report = new GATKReport();

        report.addTable(reportName, reportDescription);
        report.getTable(reportName).addPrimaryKey("cycle");

        for (SAMReadGroupRecord rg : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            String readGroupId = rg.getReadGroupId();

            report.getTable(reportName).addColumn("mismatches." + readGroupId, 0, false);
            report.getTable(reportName).addColumn(   "qualsum." + readGroupId, 0, false);
            report.getTable(reportName).addColumn(    "counts." + readGroupId, 0, false);
            report.getTable(reportName).addColumn( "errorrate." + readGroupId, 0.0f);
            report.getTable(reportName).addColumn(   "qualavg." + readGroupId, 0.0f);
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<Integer> offsets = context.getOffsets();
        List<SAMRecord> reads = context.getReads();

        for (int i = 0; i < offsets.size(); i++) {
            int offset = offsets.get(i);

            if (reads.get(i).getMappingQuality() >= MIN_MAPPING_QUAL && reads.get(i).getBaseQualities()[offset] >= MIN_BASE_QUAL) {
                char readBase = reads.get(i).getReadString().charAt(offset);

                int refIndex = ref.getBaseIndex();
                int readIndex = BaseUtils.simpleBaseToBaseIndex(readBase);

                if (!reads.get(i).getReadNegativeStrandFlag() && (!reads.get(i).getReadPairedFlag() || reads.get(i).getFirstOfPairFlag())) {
                    String readGroupId = reads.get(i).getReadGroup().getReadGroupId();

                    if (refIndex != readIndex) {
                        report.getTable(reportName).increment(offset, "mismatches." + readGroupId);
                    }
                    report.getTable(reportName).add(offset, "qualsum." + readGroupId, (int) reads.get(i).getBaseQualities()[offset]);
                    report.getTable(reportName).increment(offset, "counts." + readGroupId);
                }
            }
        }

        return null;
    }

    public Integer reduceInit() { return null; }

    public Integer reduce(Integer value, Integer sum) { return null; }

    public void onTraversalDone(Integer sum) {
        for (SAMReadGroupRecord rg : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            String readGroupId = rg.getReadGroupId();
            
            report.getTable(reportName).divideColumns("errorrate." + readGroupId, "mismatches." + readGroupId, "counts." + readGroupId);
            report.getTable(reportName).divideColumns(  "qualavg." + readGroupId,    "qualsum." + readGroupId, "counts." + readGroupId);
        }

        report.print(out);
    }
}
