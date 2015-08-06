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

package org.broadinstitute.gatk.tools.walkers.diagnostics;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportTable;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.PrintStream;

/**
 * Compute the read error rate per position
 *
 * <p>This tool computes the read error rate per position in sequence reads. It does this in the original 5'->3'
 * orientation that the read had coming off the machine. It then emits a GATKReport containing readgroup, cycle,
 * mismatches, counts, qual, and error rate for each read group in the input BAMs.</p>
 *
 * <h3>Input</h3>
 *  <p>
 *      Any number of BAM files
 *  </p>
 *
 * <h3>Output</h3>
 *  <p>
 *      A GATKReport containing readgroup, cycle, mismatches, counts, qual, and error rate.
 *
 *      For example, running this tool on the NA12878 data sets yields the following table:
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
 * <h3>Usage example</h3>
 *  <pre>
 *    java -jar GenomeAnalysisTK.jar \
 *      -T ErrorRatePerCycle \
 *      -R reference.fasta \
 *      -I my_sequence_reads.bam \
 *      -o error_rates.gatkreport.txt
 *  </pre>
 *
 * <h3>Caveat</h3>
 *
 * <p>When it is run on paired-end sequence data, this tool only uses the first read in a pair.</p>
 *
 * @author Kiran Garimella, Mark DePristo
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class ErrorRatePerCycle extends LocusWalker<Integer, Integer> {
    @Output PrintStream out;
    @Argument(fullName="min_base_quality_score", shortName="mbq", doc="Minimum base quality required to consider a base for calling", required=false)
    public Integer MIN_BASE_QUAL = 0;
    @Argument(fullName="min_mapping_quality_score", shortName="mmq", doc="Minimum read mapping quality required to consider a read for calling", required=false)
    public Integer MIN_MAPPING_QUAL = 20;

    private GATKReport report;
    private GATKReportTable table;
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
        report.addTable(reportName, reportDescription, 6, GATKReportTable.TableSortingWay.SORT_BY_ROW);
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
            final int qual = QualityUtils.errorProbToQual(errorRate);
            table.set(key, "qual", qual);
            table.set(key, "errorrate", errorRate);
        }

        report.print(out);
    }
}
