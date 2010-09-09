/*
 * Copyright (c) 2009 The Broad Institute
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use,
 *  copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following
 *  conditions:
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *  * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.UserError;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.EnumMap;
import java.util.Map;
import java.io.File;
import java.io.PrintStream;
import java.io.FileNotFoundException;


/**
 * Emits a data file containing information about callable, uncallable, poorly mapped, and other parts of the genome
 *
 * @Author depristo
 * @Date May 7, 2010
 */
@By(DataSource.REFERENCE)
public class CallableLociWalker extends LocusWalker<CallableLociWalker.CallableBaseState, CallableLociWalker.Integrator> {
    @Output
    PrintStream out;

    @Argument(fullName = "maxLowMAPQ", shortName = "mlmq", doc = "Maximum value for MAPQ to be considered a problematic mapped read.  The gap between this value and mmq are reads that are not sufficiently well mapped for calling but aren't indicative of mapping problems.", required = false)
    byte maxLowMAPQ = 1;

    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth. Defaults to 50.", required = false)
    byte minMappingQuality = 10;

    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth. Defaults to 20.", required = false)
    byte minBaseQuality = 20;

    @Argument(fullName = "minDepth", shortName = "minDepth", doc = "Minimum QC+ read depth before a locus is considered callable", required = false)
    int minDepth = 4;

    @Argument(fullName = "maxDepth", shortName = "maxDepth", doc = "Maximum read depth before a locus is considered poorly mapped", required = false)
    int maxDepth = -1;

    @Argument(fullName = "minDepthForLowMAPQ", shortName = "mdflmq", doc = "Minimum read depth before a locus is considered a potential candidate for poorly mapped", required = false)
    int minDepthLowMAPQ = 10;

    @Argument(fullName = "maxFractionOfReadsWithLowMAPQ", shortName = "frlmq", doc = "Maximum read depth before a locus is considered poorly mapped", required = false)
    double maxLowMAPQFraction = 0.1;

    @Argument(fullName = "format", shortName = "format", doc = "Output format for the system: either BED or STATE_PER_BASE", required = false)
    OutputFormat outputFormat;

    @Argument(fullName = "summary", shortName = "summary", doc = "Name of file for output summary", required = true)
    File summaryFile;

    public enum OutputFormat { BED, STATE_PER_BASE }

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD WALKER METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public boolean includeReadsWithDeletionAtLoci() { return true; }

    public void initialize() {
        try {
            PrintStream summaryOut = new PrintStream(summaryFile);
            summaryOut.close();
        } catch (FileNotFoundException e) {
            throw new UserError.CouldNotCreateOutputFile(summaryFile, e);
        }
    }

    public static class Integrator {
        long counts[] = new long[CalledState.values().length];
        CallableBaseState state = null;
    }

    public static class CallableBaseState {
        public GenomeLoc loc;
        public CalledState state;

        public CallableBaseState(GenomeLoc loc, CalledState state) {
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
        public boolean changingState( CalledState newState ) {
            return state != newState;
        }

        /**
         * Updating the location of this CalledBaseState by the new stop location
         * @param newStop
         */
        public void update(GenomeLoc newStop) {
            loc = GenomeLocParser.createGenomeLoc(loc.getContigIndex(), loc.getStart(), newStop.getStop());
        }

        public String toString() {
            return String.format("%s %d %d %s", loc.getContig(), loc.getStart(), loc.getStop(), state);
            //return String.format("%s %d %d %d %s", loc.getContig(), loc.getStart(), loc.getStop(), loc.getStop() - loc.getStart() + 1, state);
        }
    }

    public enum CalledState { REF_N, CALLABLE, NO_COVERAGE, LOW_COVERAGE, EXCESSIVE_COVERAGE, POOR_MAPPING_QUALITY }

    public Integrator reduceInit() {
        return new Integrator();
    }

    public CallableBaseState map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        CalledState state;

        if ( BaseUtils.isNBase(ref.getBase()) ) {
            state = CalledState.REF_N;
        } else {
            // count up the depths of all and QC+ bases
            int rawDepth = 0, QCDepth = 0, lowMAPQDepth = 0;
            for (PileupElement e : context.getBasePileup()) {
                rawDepth++;

                if ( e.getMappingQual() <= maxLowMAPQ )
                    lowMAPQDepth++;

                if ( e.getMappingQual() >= minMappingQuality && ( e.getQual() >= minBaseQuality || e.isDeletion() ) ) {
                    QCDepth++;
                }
            }

            //System.out.printf("%s rawdepth = %d QCDepth = %d lowMAPQ = %d%n", context.getLocation(), rawDepth, QCDepth, lowMAPQDepth);
            if ( rawDepth == 0 ) {
                state = CalledState.NO_COVERAGE;
            } else if ( rawDepth >= minDepthLowMAPQ && MathUtils.ratio( lowMAPQDepth, rawDepth ) >= maxLowMAPQFraction ) {
                state = CalledState.POOR_MAPPING_QUALITY;
            } else if ( QCDepth < minDepth ) {
                state = CalledState.LOW_COVERAGE;
            } else if ( rawDepth >= maxDepth && maxDepth != -1 ) {
                state = CalledState.EXCESSIVE_COVERAGE;
            } else {
                state = CalledState.CALLABLE;
            }
        }

        return new CallableBaseState(context.getLocation(), state);
    }

    public Integrator reduce(CallableBaseState state, Integrator integrator) {
        // update counts
        integrator.counts[state.getState().ordinal()]++;

        if ( outputFormat == OutputFormat.STATE_PER_BASE ) {
            out.println(state.toString());
        } else {
            // format is integrating
            if ( integrator.state == null )
                integrator.state = state;
            else if ( state.getLocation().getStart() != integrator.state.getLocation().getStop() + 1 ||
                    integrator.state.changingState(state.getState()) ) {
                out.println(integrator.state.toString());
                integrator.state = state;
            } else {
                integrator.state.update(state.getLocation());
            }
        }

        return integrator;
    }


    ////////////////////////////////////////////////////////////////////////////////////
    // INTERVAL ON TRAVERSAL DONE
    ////////////////////////////////////////////////////////////////////////////////////

    public void onTraversalDone(Integrator result) {
        // print out the last state
        if ( result != null ) {
            out.println(result.state.toString());

            try {
                PrintStream summaryOut = new PrintStream(summaryFile);
                summaryOut.printf("%30s %s%n", "state", "nBases");
                for ( CalledState state : CalledState.values() ) {
                    summaryOut.printf("%30s %d%n", state, result.counts[state.ordinal()]);
                }

                summaryOut.close();
            } catch (FileNotFoundException e) {
                throw new UserError.CouldNotCreateOutputFile(summaryFile, e);
            }
        }
    }
}