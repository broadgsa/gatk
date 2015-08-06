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

package org.broadinstitute.gatk.tools.walkers.coverage;

import org.broadinstitute.gatk.utils.commandline.Advanced;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.PileupElement;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;


/**
 * Collect statistics on callable, uncallable, poorly mapped, and other parts of the genome
 *
 * <p>
 * A very common question about a NGS set of reads is what areas of the genome are considered callable. This tool
 * considers the coverage at each locus and emits either a per base state or a summary interval BED file that
 * partitions the genomic intervals into the following callable states:
 * <dl>
 * <dt>REF_N</dt>
 * <dd>The reference base was an N, which is not considered callable the GATK</dd>
 * <dt>PASS</dt>
 * <dd>The base satisfied the min. depth for calling but had less than maxDepth to avoid having EXCESSIVE_COVERAGE</dd>
 * <dt>NO_COVERAGE</dt>
 * <dd>Absolutely no reads were seen at this locus, regardless of the filtering parameters</dd>
 * <dt>LOW_COVERAGE</dt>
 * <dd>There were fewer than min. depth bases at the locus, after applying filters</dd>
 * <dt>EXCESSIVE_COVERAGE</dt>
 * <dd>More than -maxDepth read at the locus, indicating some sort of mapping problem</dd>
 * <dt>POOR_MAPPING_QUALITY</dt>
 * <dd>More than --maxFractionOfReadsWithLowMAPQ at the locus, indicating a poor mapping quality of the reads</dd>
 * </dl>
 * </p>
 * <p/>
 * <h3>Input</h3>
 * <p>
 * A BAM file containing <b>exactly one sample</b>.
 * </p>
 * <p/>
 * <h3>Output</h3>
 * <p>
 *     A file with the callable status covering each base and a table of callable status x count of all examined bases
 * </p>
 * <h3>Usage example</h3>
 * <pre>
 *  java -jar GenomeAnalysisTK.jar \
 *     -T CallableLoci \
 *     -R reference.fasta \
 *     -I myreads.bam \
 *     -summary table.txt \
 *     -o callable_status.bed
 * </pre>
 * <p/>
 * would produce a BED file that looks like:
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
 * </pre>
 * as well as a summary table that looks like:
 * <p/>
 * <pre>
 *                        state nBases
 *                        REF_N 0
 *                         PASS 996046
 *                  NO_COVERAGE 121
 *                 LOW_COVERAGE 928
 *           EXCESSIVE_COVERAGE 0
 *         POOR_MAPPING_QUALITY 2906
 * </pre>
 *
 * @author Mark DePristo
 * @since May 7, 2010
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@By(DataSource.REFERENCE)
public class CallableLoci extends LocusWalker<CallableLoci.CallableBaseState, CallableLoci.Integrator> {
    @Output
    PrintStream out;

    /**
     * Callable loci summary counts will be written to this file.
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
     * The output of this tool will be written in this format.  The recommended option is BED.
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
         * Emit chr start stop state quads for each base.  Produces a potentially disastrously
         * large amount of output.
         */
        STATE_PER_BASE
    }

    public enum CalledState {
        /**
         * the reference base was an N, which is not considered callable by the GATK
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
         * there were fewer than min. depth bases at the locus, after applying filters
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