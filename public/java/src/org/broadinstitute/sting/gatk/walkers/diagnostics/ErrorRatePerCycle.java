package org.broadinstitute.sting.gatk.walkers.diagnostics;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTableV2;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;

/**
 * Computes the read error rate per position in read (in the original 5'->3' orientation that the read had coming off the machine)
 *
 * Emits a GATKReport containing readgroup, cycle, mismatches, counts, qual, and error rate for each read
 * group in the input BAMs FOR ONLY THE FIRST OF PAIR READS.
 *
 * <h2>Input</h2>
 *  <p>
 *      Any number of BAM files
 *  </p>
 *
 * <h2>Output</h2>
 *  <p>
 *      GATKReport containing readgroup, cycle, mismatches, counts, qual, and error rate.
 *
 *      For example, running this tool on the NA12878 data sets:
 *
 *      <pre>
 *      ##:GATKReport.v0.2 ErrorRatePerCycle : The error rate per sequenced position in the reads
 *      readgroup  cycle  mismatches  counts  qual  errorrate
 *      20FUK.1        0          80   23368    25   3.47e-03
 *      20FUK.1        1          40   23433    28   1.75e-03
 *      20FUK.1        2          36   23453    28   1.58e-03
 *      20FUK.1        3          26   23476    29   1.15e-03
 *      20FUK.1        4          32   23495    29   1.40e-03
 *      up to 101 cycles
 *      20FUK.2        0          77   20886    24   3.73e-03
 *      20FUK.2        1          28   20920    29   1.39e-03
 *      20FUK.2        2          24   20931    29   1.19e-03
 *      20FUK.2        3          30   20940    28   1.48e-03
 *      20FUK.2        4          25   20948    29   1.24e-03
 *      up to 101 cycles
 *      20FUK.3        0          78   22038    24   3.58e-03
 *      20FUK.3        1          40   22091    27   1.86e-03
 *      20FUK.3        2          23   22108    30   1.09e-03
 *      20FUK.3        3          36   22126    28   1.67e-03
 *      </pre>
 *  </p>
 *
 * <h2>Examples</h2>
 *  <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T ErrorRatePerCycle
 *      -I bundle/current/b37/NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam
 *      -R bundle/current/b37/human_g1k_v37.fasta
 *      -o example.gatkreport.txt
 *  </pre>
 *
 * @author Kiran Garimella, Mark DePristo
 */
public class ErrorRatePerCycle extends LocusWalker<Integer, Integer> {
    @Output PrintStream out;
    @Argument(fullName="min_base_quality_score", shortName="mbq", doc="Minimum base quality required to consider a base for calling", required=false)
    public Integer MIN_BASE_QUAL = 0;
    @Argument(fullName="min_mapping_quality_score", shortName="mmq", doc="Minimum read mapping quality required to consider a read for calling", required=false)
    public Integer MIN_MAPPING_QUAL = 20;

    private GATKReport report;
    private GATKReportTableV2 table;
    private final static String reportName = "ErrorRatePerCycle";
    private final static String reportDescription = "The error rate per sequenced position in the reads";

    /**
     * Allows us to use multiple records for the key (read group x cycle)
     */
    private static class TableKey implements Comparable<TableKey> {
        final String readGroup;
        final int cycle;

        private TableKey(final String readGroup, final int cycle) {
            this.readGroup = readGroup;
            this.cycle = cycle;
        }

        // Must overload hashCode and equals to properly work with GATKReportColumn
        @Override
        public int hashCode() {
            return readGroup.hashCode() + 33 * cycle;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final TableKey oKey = (TableKey) o;

            if ( cycle != oKey.cycle ) return false;
            if ( !readGroup.equals(oKey.readGroup) ) return false;

            return true;
        }

        @Override
        public int compareTo(final TableKey tableKey) {
            final int scmp = readGroup.compareTo(tableKey.readGroup);
            if ( scmp == 0 )
                return Integer.valueOf(cycle).compareTo(tableKey.cycle);
            else
                return scmp;
        }
    }

    public void initialize() {
        report = new GATKReport();
        report.addTable(reportName, reportDescription, 6, true);
        table = report.getTable(reportName);
        table.addColumn("readgroup");
        table.addColumn("cycle");
        table.addColumn("mismatches");
        table.addColumn("counts");
        table.addColumn("qual");
        table.addColumn("errorrate", "%.2e");
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        for ( final PileupElement p : context.getBasePileup() ) {
            final GATKSAMRecord read = p.getRead();
            final int offset = p.getOffset();
            final boolean firstOfPair = ! read.getReadPairedFlag() || read.getFirstOfPairFlag();

            if ( firstOfPair && read.getMappingQuality() >= MIN_MAPPING_QUAL && p.getQual() >= MIN_BASE_QUAL ) {
                final byte readBase = p.getBase();
                final byte refBase = ref.getBase();
                final int cycle = offset;

                if ( BaseUtils.isRegularBase(readBase) && BaseUtils.isRegularBase(refBase) ) {
                    final TableKey key = new TableKey(read.getReadGroup().getReadGroupId(), cycle);

                    if ( ! table.containsRowID(key) ) {
                        table.set(key, "cycle", cycle);
                        table.set(key, "readgroup", read.getReadGroup().getReadGroupId());
                        table.set(key, "counts", 0);
                        table.set(key, "mismatches", 0);
                    }

                    table.increment(key, "counts");
                    if (readBase != refBase)
                        table.increment(key, "mismatches");
                }
            }
        }

        return null;
    }

    public Integer reduceInit() { return null; }

    public Integer reduce(Integer value, Integer sum) { return null; }

    public void onTraversalDone(Integer sum) {
        for ( Object key : table.getRowIDs() ) {
            final int mismatches = (Integer)table.get(key, "mismatches");
            final int count = (Integer)table.get(key, "counts");
            final double errorRate = (mismatches + 1) / (1.0*(count + 1));
            final int qual = QualityUtils.probToQual(1-errorRate, 0.0);
            table.set(key, "qual", qual);
            table.set(key, "errorrate", errorRate);
        }

        report.print(out);
    }
}
