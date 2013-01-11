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

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.Median;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.text.DateFormat;
import java.util.HashMap;
import java.util.Map;

/**
 * Emits a GATKReport containing read group, sample, library, platform, center, sequencing data,
 * paired end status, simple read type name (e.g. 2x76) median insert size and median read length
 * for each read group in every provided BAM file
 *
 * Note that this walker stops when all read groups have been observed at least a few thousand times so that
 * the median statistics are well determined.  It is safe to run it WG and it'll finish in an appropriate
 * timeframe.
 *
 * <h2>Input</h2>
 *  <p>
 *      Any number of BAM files
 *  </p>
 *
 * <h2>Output</h2>
 *  <p>
 *      GATKReport containing read group, sample, library, platform, center, median insert size and median read length.
 *
 *      For example, running this tool on the NA12878 data sets:
 *
 *      <pre>
 *      ##:GATKReport.v0.2 ReadGroupProperties : Table of read group properties
 *      readgroup  sample   library       platform  center  date     has.any.reads  is.paired.end  n.reads.analyzed  simple.read.type  median.read.length  median.insert.size
 *      20FUK.1    NA12878  Solexa-18483  illumina  BI      2/2/10   true           true                        498  2x101                            101                 386
 *      20FUK.2    NA12878  Solexa-18484  illumina  BI      2/2/10   true           true                        476  2x101                            101                 417
 *      20FUK.3    NA12878  Solexa-18483  illumina  BI      2/2/10   true           true                        407  2x101                            101                 387
 *      20FUK.4    NA12878  Solexa-18484  illumina  BI      2/2/10   true           true                        389  2x101                            101                 415
 *      20FUK.5    NA12878  Solexa-18483  illumina  BI      2/2/10   true           true                        433  2x101                            101                 386
 *      20FUK.6    NA12878  Solexa-18484  illumina  BI      2/2/10   true           true                        480  2x101                            101                 418
 *      20FUK.7    NA12878  Solexa-18483  illumina  BI      2/2/10   true           true                        450  2x101                            101                 386
 *      20FUK.8    NA12878  Solexa-18484  illumina  BI      2/2/10   true           true                        438  2x101                            101                 418
 *      20GAV.1    NA12878  Solexa-18483  illumina  BI      1/26/10  true           true                        490  2x101                            101                 391
 *      20GAV.2    NA12878  Solexa-18484  illumina  BI      1/26/10  true           true                        485  2x101                            101                 417
 *      20GAV.3    NA12878  Solexa-18483  illumina  BI      1/26/10  true           true                        460  2x101                            101                 392
 *      20GAV.4    NA12878  Solexa-18484  illumina  BI      1/26/10  true           true                        434  2x101                            101                 415
 *      20GAV.5    NA12878  Solexa-18483  illumina  BI      1/26/10  true           true                        479  2x101                            101                 389
 *      20GAV.6    NA12878  Solexa-18484  illumina  BI      1/26/10  true           true                        461  2x101                            101                 416
 *      20GAV.7    NA12878  Solexa-18483  illumina  BI      1/26/10  true           true                        509  2x101                            101                 386
 *      20GAV.8    NA12878  Solexa-18484  illumina  BI      1/26/10  true           true                        476  2x101                            101                 410                           101                 414
 *      </pre>
 *  </p>
 *
 * <h2>Examples</h2>
 *  <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T ReadGroupProperties
 *      -I example1.bam -I example2.bam etc
 *      -R reference.fasta
 *      -o example.gatkreport.txt
 *  </pre>
 *
 * @author Mark DePristo
 */
@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
public class ReadGroupProperties extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Argument(shortName="maxElementsForMedian", doc="Calculate median from the first maxElementsForMedian values observed", required=false)
    public int MAX_VALUES_FOR_MEDIAN = 10000;

    private final static String TABLE_NAME = "ReadGroupProperties";
    private final Map<String, PerReadGroupInfo> readGroupInfo = new HashMap<String, PerReadGroupInfo>();

    private class PerReadGroupInfo {
        public final Median<Integer> readLength = new Median<Integer>(MAX_VALUES_FOR_MEDIAN);
        public final Median<Integer> insertSize = new Median<Integer>(MAX_VALUES_FOR_MEDIAN);
        public int nReadsSeen = 0, nReadsPaired = 0;

        public boolean needsMoreData() {
            return ! readLength.isFull() || ! insertSize.isFull();
        }
    }

    @Override
    public void initialize() {
        for ( final SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
            readGroupInfo.put(rg.getId(), new PerReadGroupInfo());
        }
    }

    @Override
    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        return ! (read.getReadFailsVendorQualityCheckFlag() || read.getReadUnmappedFlag());
    }

    @Override
    public boolean isDone() {
        for ( PerReadGroupInfo info : readGroupInfo.values() ) {
            if ( info.needsMoreData() )
                return false;
        }

        return true;
    }

    @Override
    public Integer map(ReferenceContext referenceContext, GATKSAMRecord read, RefMetaDataTracker RefMetaDataTracker) {
        final String rgID = read.getReadGroup().getId();
        final PerReadGroupInfo info = readGroupInfo.get(rgID);

        if ( info.needsMoreData() ) {
            info.readLength.add(read.getReadLength());
            info.nReadsSeen++;
            if ( read.getReadPairedFlag() ) {
                info.nReadsPaired++;
                if ( read.getInferredInsertSize() != 0) {
                    info.insertSize.add(Math.abs(read.getInferredInsertSize()));
                }
            }
        }

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer integer, Integer integer1) {
        return null;
    }

    @Override
    public void onTraversalDone(Integer sum) {
        final GATKReport report = new GATKReport();
        report.addTable(TABLE_NAME, "Table of read group properties", 12);
        GATKReportTable table = report.getTable(TABLE_NAME);
        DateFormat dateFormatter = DateFormat.getDateInstance(DateFormat.SHORT);

        table.addColumn("readgroup");
        //* Emits a GATKReport containing read group, sample, library, platform, center, median insert size and
        //* median read length for each read group in every BAM file.
        table.addColumn("sample", "%s");
        table.addColumn("library", "%s");
        table.addColumn("platform", "%s");
        table.addColumn("center", "%s");
        table.addColumn("date", "%s");
        table.addColumn("has.any.reads");
        table.addColumn("is.paired.end");
        table.addColumn("n.reads.analyzed", "%d");
        table.addColumn("simple.read.type", "%s");
        table.addColumn("median.read.length");
        table.addColumn("median.insert.size");

        for ( final SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
            final String rgID = rg.getId();
            table.addRowID(rgID, true);
            PerReadGroupInfo info = readGroupInfo.get(rgID);

            // we are paired if > 25% of reads are paired
            final boolean isPaired = info.nReadsPaired / (1.0 * (info.nReadsSeen+1)) > 0.25;
            final boolean hasAnyReads = info.nReadsSeen > 0;
            final int readLength = info.readLength.getMedian(0);

            setTableValue(table, rgID, "sample", rg.getSample());
            setTableValue(table, rgID, "library", rg.getLibrary());
            setTableValue(table, rgID, "platform", rg.getPlatform());
            setTableValue(table, rgID, "center", rg.getSequencingCenter());
            try {
                setTableValue(table, rgID, "date", rg.getRunDate() != null ? dateFormatter.format(rg.getRunDate()) : "NA");
            } catch ( NullPointerException e ) {
                // TODO: remove me when bug in Picard is fixed that causes NPE when date isn't present
                setTableValue(table, rgID, "date", "NA");
            }
            setTableValue(table, rgID, "has.any.reads", hasAnyReads);
            setTableValue(table, rgID, "is.paired.end", isPaired);
            setTableValue(table, rgID, "n.reads.analyzed", info.nReadsSeen);
            setTableValue(table, rgID, "simple.read.type", hasAnyReads ? String.format("%dx%d", isPaired ? 2 : 1, readLength) : "NA");
            setTableValue(table, rgID, "median.read.length", hasAnyReads ? readLength : "NA" );
            setTableValue(table, rgID, "median.insert.size", hasAnyReads && isPaired ? info.insertSize.getMedian(0) : "NA" );
        }

        report.print(out);
    }

    private final void setTableValue(GATKReportTable table, final String rgID, final String key, final Object value) {
        table.set(rgID, key, value == null ? "NA" : value);
    }
}
