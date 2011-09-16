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

import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * Test routine for new VariantContext object
 */
public class CompareCallableLociWalker extends RodWalker<List<CallableLociWalker.CallableBaseState>, long[][]> {
    @Output
    protected PrintStream out;

    @Input(fullName="comp1", shortName = "comp1", doc="First comparison track name", required=true)
    public RodBinding<BEDFeature> compTrack1;

    @Input(fullName="comp2", shortName = "comp2", doc="Second comparison track name", required=true)
    public RodBinding<BEDFeature> compTrack2;

    @Argument(shortName="printState", doc="If provided, prints sites satisfying this state pair", required=false)
    protected String printState = null;

    CallableLociWalker.CalledState printState1 = CallableLociWalker.CalledState.REF_N;
    CallableLociWalker.CalledState printState2 = CallableLociWalker.CalledState.REF_N;

    // --------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    // --------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if ( printState != null ) {
            String[] states = printState.split(",");
            printState1 = CallableLociWalker.CalledState.valueOf(states[0]);
            printState2 = CallableLociWalker.CalledState.valueOf(states[1]);
        }
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    // --------------------------------------------------------------------------------------------------------------
    public List<CallableLociWalker.CallableBaseState> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker != null ) {
            CallableLociWalker.CallableBaseState comp1 = getCallableBaseState(tracker, compTrack1);
            CallableLociWalker.CallableBaseState comp2 = getCallableBaseState(tracker, compTrack2);

            if ( printState != null && comp1.getState() == printState1 && comp2.getState() == printState2 ) {
                out.printf("%s %s %s %s%n", comp1.getLocation(), comp1.getState(), comp2.getLocation(), comp2.getState());
            }

            return Arrays.asList(comp1, comp2);
        } else {
            return null;
        }
    }

    private CallableLociWalker.CallableBaseState getCallableBaseState(RefMetaDataTracker tracker, RodBinding<BEDFeature> rodBinding) {
        //System.out.printf("tracker %s%n", tracker);
        List<BEDFeature> bindings = tracker.getValues(rodBinding);
        if ( bindings.size() != 1 ) {
            throw new UserException.MalformedFile(String.format("%s track isn't a properly formated CallableBases object!", rodBinding.getName()));
        }

        BEDFeature bed = bindings.get(0);
        GenomeLoc loc = getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), bed.getStart(), bed.getEnd());
        CallableLociWalker.CalledState state = CallableLociWalker.CalledState.valueOf(bed.getName());
        return new CallableLociWalker.CallableBaseState(getToolkit().getGenomeLocParser(),loc, state);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    // --------------------------------------------------------------------------------------------------------------
    public long[][] reduceInit() {
        int n = CallableLociWalker.CalledState.values().length;
        return new long[n][n];
    }

    public long[][] reduce(List<CallableLociWalker.CallableBaseState> comps, long[][] sum) {
        if ( comps != null ) {
            CallableLociWalker.CallableBaseState comp1 = comps.get(0);
            CallableLociWalker.CallableBaseState comp2 = comps.get(1);

            sum[comp1.getState().ordinal()][comp2.getState().ordinal()]++;
        }

        return sum;
    }

    public void onTraversalDone(long[][] result) {
        for ( CallableLociWalker.CalledState state1 : CallableLociWalker.CalledState.values() ) {
            for ( CallableLociWalker.CalledState state2 : CallableLociWalker.CalledState.values() ) {
                out.printf("%s %s %s %s %d%n", compTrack1.getName(), compTrack2.getName(), state1, state2, result[state1.ordinal()][state2.ordinal()]);
            }
        }
    }
}