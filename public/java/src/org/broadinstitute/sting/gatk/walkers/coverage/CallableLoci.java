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

package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.variant.utils.BaseUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;


/**
 * Emits a data file containing information about callable, uncallable, poorly mapped, and other parts of the genome
 * <p/>
 * <p>
 * A very common question about a NGS set of reads is what areas of the genome are considered callable. The system
 * considers the coverage at each locus and emits either a per base state or a summary interval BED file that
 * partitions the genomic intervals into the following callable states:
 * <dl>
 * <dt>REF_N</dt>
 * <dd>the reference base was an N, which is not considered callable the GATK</dd>
 * <dt>PASS</dt>
 * <dd>the base satisfied the min. depth for calling but had less than maxDepth to avoid having EXCESSIVE_COVERAGE</dd>
 * <dt>NO_COVERAGE</dt>
 * <dd>absolutely no reads were seen at this locus, regardless of the filtering parameters</dd>
 * <dt>LOW_COVERAGE</dt>
 * <dd>there were less than min. depth bases at the locus, after applying filters</dd>
 * <dt>EXCESSIVE_COVERAGE</dt>
 * <dd>more than -maxDepth read at the locus, indicating some sort of mapping problem</dd>
 * <dt>POOR_MAPPING_QUALITY</dt>
 * <dd>more than --maxFractionOfReadsWithLowMAPQ at the locus, indicating a poor mapping quality of the reads</dd>
 * </dl>
 * </p>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * A BAM file containing <b>exactly one sample</b>.
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * <ul>
 * <li>-o: a OutputFormatted (recommended BED) file with the callable status covering each base</li>
 * <li>-summary: a table of callable status x count of all examined bases</li>
 * </ul>
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *     -T CallableLociWalker \
 *     -I my.bam \
 *     -summary my.summary \
 *     -o my.bed
 * </pre>
 * <p/>
 * would produce a BED file (my.bed) that looks like:
 * <p/>
 * <pre>
 *     20 10000000 10000864 PASS
 *     20 10000865 10000985 POOR_MAPPING_QUALITY
 *     20 10000986 10001138 PASS
 *     20 10001139 10001254 POOR_MAPPING_QUALITY
 *     20 10001255 10012255 PASS
 *     20 10012256 10012259 POOR_MAPPING_QUALITY
 *     20 10012260 10012263 PASS
 *     20 10012264 10012328 POOR_MAPPING_QUALITY
 *     20 10012329 10012550 PASS
 *     20 10012551 10012551 LOW_COVERAGE
 *     20 10012552 10012554 PASS
 *     20 10012555 10012557 LOW_COVERAGE
 *     20 10012558 10012558 PASS
 *     et cetera...
 * </pre>
 * as well as a summary table that looks like:
 * <p/>
 * <pre>
 *                        state nBases
 *                        REF_N 0
 *                     PASS 996046
 *                  NO_COVERAGE 121
 *                 LOW_COVERAGE 928
 *           EXCESSIVE_COVERAGE 0
 *         POOR_MAPPING_QUALITY 2906
 * </pre>
 *
 * @author Mark DePristo
 * @since May 7, 2010
 */
@DocumentedGATKFeature( groupName = "BAM Processing and Analysis Tools", extraDocs = {CommandLineGATK.class} )
@By(DataSource.REFERENCE)
public class CallableLoci extends LocusWalker<CallableLoci.CallableBaseState, CallableLoci.Integrator> {
    @Output
    PrintStream out;

    /**
     * Callable loci summary counts (see outputs) will be written to this file.
     */
    @Output(fullName = "summary", shortName = "summary", doc = "Name of file for output summary", required = true)
    File summaryFile;

    /**
     * The gap between this value and mmq are reads that are not sufficiently well mapped for calling but
     * aren't indicative of mapping problems.  For example, if maxLowMAPQ = 1 and mmq = 20, then reads with
     * MAPQ == 0 are poorly mapped, MAPQ >= 20 are considered as contributing to calling, where
     * reads with MAPQ >= 1 and < 20 are not bad in and of themselves but aren't sufficiently good to contribute to
     * calling.  In effect this reads are invisible, driving the base to the NO_ or LOW_COVERAGE states
     */
    @Argument(fullName = "maxLowMAPQ", shortName = "mlmq", doc = "Maximum value for MAPQ to be considered a problematic mapped read.", required = false)
    byte maxLowMAPQ = 1;

    /**
     * Reads with MAPQ > minMappingQuality are treated as usable for variation detection, contributing to the PASS
     * state.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth.", required = false)
    byte minMappingQuality = 10;

    /**
     * Bases with less than minBaseQuality are viewed as not sufficiently high quality to contribute to the PASS state
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth.", required = false)
    byte minBaseQuality = 20;

    /**
     * If the number of QC+ bases (on reads with MAPQ > minMappingQuality and with base quality > minBaseQuality) exceeds this
     * value and is less than maxDepth the site is considered PASS.
     */
    @Advanced
    @Argument(fullName = "minDepth", shortName = "minDepth", doc = "Minimum QC+ read depth before a locus is considered callable", required = false)
    int minDepth = 4;

    /**
     * If the QC+ depth exceeds this value the site is considered to have EXCESSIVE_DEPTH
     */
    @Argument(fullName = "maxDepth", shortName = "maxDepth", doc = "Maximum read depth before a locus is considered poorly mapped", required = false)
    int maxDepth = -1;

    /**
     * We don't want to consider a site as POOR_MAPPING_QUALITY just because it has two reads, and one is MAPQ.  We
     * won't assign a site to the POOR_MAPPING_QUALITY state unless there are at least minDepthForLowMAPQ reads
     * covering the site.
     */
    @Advanced
    @Argument(fullName = "minDepthForLowMAPQ", shortName = "mdflmq", doc = "Minimum read depth before a locus is considered a potential candidate for poorly mapped", required = false)
    int minDepthLowMAPQ = 10;

    /**
     * If the number of reads at this site is greater than minDepthForLowMAPQ and the fraction of reads with low mapping quality
     * exceeds this fraction then the site has POOR_MAPPING_QUALITY.
     */
    @Argument(fullName = "maxFractionOfReadsWithLowMAPQ", shortName = "frlmq", doc = "If the fraction of reads at a base with low mapping quality exceeds this value, the site may be poorly mapped", required = false)
    double maxLowMAPQFraction = 0.1;

    /**
     * The output of this walker will be written in this format.  The recommended option is BED.
     */
    @Advanced
    @Argument(fullName = "format", shortName = "format", doc = "Output format", required = false)
    OutputFormat outputFormat = OutputFormat.BED;

    public enum OutputFormat {
        /**
         * The output will be written as a BED file.  There's a BED element for each
         * continuous run of callable states (i.e., PASS, REF_N, etc).  This is the recommended
         * format
         */
        BED,

        /**
         * Emit chr start stop state quads for each base.  Produces a potentially disasterously
         * large amount of output.
         */
        STATE_PER_BASE
    }

    public enum CalledState {
        /**
         * the reference base was an N, which is not considered callable the GATK
         */
        REF_N,
        /**
         * the base satisfied the min. depth for calling but had less than maxDepth to avoid having EXCESSIVE_COVERAGE
         */
        CALLABLE,
        /**
         * absolutely no reads were seen at this locus, regardless of the filtering parameters
         */
        NO_COVERAGE,
        /**
         * there were less than min. depth bases at the locus, after applying filters
         */
        LOW_COVERAGE,
        /**
         * more than -maxDepth read at the locus, indicating some sort of mapping problem
         */
        EXCESSIVE_COVERAGE,
        /**
         * more than --maxFractionOfReadsWithLowMAPQ at the locus, indicating a poor mapping quality of the reads
         */
        POOR_MAPPING_QUALITY
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD WALKER METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    @Override
    public boolean includeReadsWithDeletionAtLoci() {
        return true;
    }

    @Override
    public void initialize() {
        if (getSampleDB().getSamples().size() != 1) {
            throw new UserException.BadArgumentValue("-I", "CallableLoci only works for a single sample, but multiple samples were found in the provided BAM files: " + getSampleDB().getSamples());
        }

        try {
            PrintStream summaryOut = new PrintStream(summaryFile);
            summaryOut.close();
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(summaryFile, e);
        }
    }

    protected static class Integrator {
        final long counts[] = new long[CalledState.values().length];
        CallableBaseState state = null;
    }

    protected static class CallableBaseState implements HasGenomeLocation {
        final public GenomeLocParser genomeLocParser;
        public GenomeLoc loc;
        final public CalledState state;

        public CallableBaseState(GenomeLocParser genomeLocParser, GenomeLoc loc, CalledState state) {
            this.genomeLocParser = genomeLocParser;
            this.loc = loc;
            this.state = state;
        }

        public GenomeLoc getLocation() {
            return loc;
        }

        public CalledState getState() {
            return state;
        }

        // update routines
        public boolean changingState(CalledState newState) {
            return state != newState;
        }

        /**
         * Updating the location of this CalledBaseState by the new stop location
         *
         * @param newStop
         */
        public void update(GenomeLoc newStop) {
            loc = genomeLocParser.createGenomeLoc(loc.getContig(), loc.getStart(), newStop.getStop());
        }

        public String toString() {
            return String.format("%s\t%d\t%d\t%s", loc.getContig(), loc.getStart()-1, loc.getStop(), state);
        }
    }

    @Override
    public CallableBaseState map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        CalledState state;

        if ( BaseUtils.isNBase(ref.getBase())) {
            state = CalledState.REF_N;
        } else {
            // count up the depths of all and QC+ bases
            int rawDepth = 0, QCDepth = 0, lowMAPQDepth = 0;
            for (PileupElement e : context.getBasePileup()) {
                rawDepth++;

                if (e.getMappingQual() <= maxLowMAPQ)
                    lowMAPQDepth++;

                if (e.getMappingQual() >= minMappingQuality && (e.getQual() >= minBaseQuality || e.isDeletion())) {
                    QCDepth++;
                }
            }

            //System.out.printf("%s rawdepth = %d QCDepth = %d lowMAPQ = %d%n", context.getLocation(), rawDepth, QCDepth, lowMAPQDepth);
            if (rawDepth == 0) {
                state = CalledState.NO_COVERAGE;
            } else if (rawDepth >= minDepthLowMAPQ && MathUtils.ratio(lowMAPQDepth, rawDepth) >= maxLowMAPQFraction) {
                state = CalledState.POOR_MAPPING_QUALITY;
            } else if (QCDepth < minDepth) {
                state = CalledState.LOW_COVERAGE;
            } else if (rawDepth >= maxDepth && maxDepth != -1) {
                state = CalledState.EXCESSIVE_COVERAGE;
            } else {
                state = CalledState.CALLABLE;
            }
        }

        return new CallableBaseState(getToolkit().getGenomeLocParser(), context.getLocation(), state);
    }

    @Override
    public Integrator reduceInit() {
        return new Integrator();
    }

    @Override
    public Integrator reduce(CallableBaseState state, Integrator integrator) {
        // update counts
        integrator.counts[state.getState().ordinal()]++;

        if (outputFormat == OutputFormat.STATE_PER_BASE) {
            out.println(state.toString());
        }

        // format is integrating
        if (integrator.state == null)
            integrator.state = state;
        else if (state.getLocation().getStart() != integrator.state.getLocation().getStop() + 1 ||
                integrator.state.changingState(state.getState())) {
            out.println(integrator.state.toString());
            integrator.state = state;
        } else {
            integrator.state.update(state.getLocation());
        }

        return integrator;
    }


    ////////////////////////////////////////////////////////////////////////////////////
    // INTERVAL ON TRAVERSAL DONE
    ////////////////////////////////////////////////////////////////////////////////////

    @Override
    public void onTraversalDone(Integrator result) {
        // print out the last state
        if (result != null) {
            if (outputFormat == OutputFormat.BED)  // get the last interval
                out.println(result.state.toString());

            try {
                PrintStream summaryOut = new PrintStream(summaryFile);
                summaryOut.printf("%30s %s%n", "state", "nBases");
                for (CalledState state : CalledState.values()) {
                    summaryOut.printf("%30s %d%n", state, result.counts[state.ordinal()]);
                }
                summaryOut.close();
            } catch (FileNotFoundException e) {
                throw new UserException.CouldNotCreateOutputFile(summaryFile, e);
            }
        }
    }
}