/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.diagnostics;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
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
@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
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
            final int qual = QualityUtils.probToQual(1-errorRate, 0.0);
            table.set(key, "qual", qual);
            table.set(key, "errorrate", errorRate);
        }

        report.print(out);
    }
}
