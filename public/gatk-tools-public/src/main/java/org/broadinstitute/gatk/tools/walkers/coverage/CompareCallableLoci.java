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

import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * Compare callability statistics
 *
 * <p>This tool can be used to evaluate how different sequence datasets compare in terms of "callability"
 * based on the output of the CallableLoci tool. </p>
 *
 *
 * <h3>Input</h3>
 * <p>
 * Two files to compare, output by two runs of CallableLoci
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A table showing the callability status of each interval of interest in the two comparison sets and whether they match.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R reference.fasta \
 *   -T CompareCallableLoci \
 *   -comp1 callable_loci_1.bed \
 *   -comp2 callable_loci_2.bed \
 *   [-L input.intervals \]
 *   -o comparison.table
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class CompareCallableLoci extends RodWalker<List<CallableLoci.CallableBaseState>, long[][]> {
    @Output
    protected PrintStream out;

    @Input(fullName="comp1", shortName = "comp1", doc="First comparison track name", required=true)
    public RodBinding<BEDFeature> compTrack1;

    @Input(fullName="comp2", shortName = "comp2", doc="Second comparison track name", required=true)
    public RodBinding<BEDFeature> compTrack2;

    @Argument(shortName="printState", doc="If provided, prints sites satisfying this state pair", required=false)
    protected String printState = null;

    CallableLoci.CalledState printState1 = CallableLoci.CalledState.REF_N;
    CallableLoci.CalledState printState2 = CallableLoci.CalledState.REF_N;

    // --------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    // --------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if ( printState != null ) {
            String[] states = printState.split(",");
            printState1 = CallableLoci.CalledState.valueOf(states[0]);
            printState2 = CallableLoci.CalledState.valueOf(states[1]);
        }
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    // --------------------------------------------------------------------------------------------------------------
    public List<CallableLoci.CallableBaseState> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker != null ) {
            CallableLoci.CallableBaseState comp1 = getCallableBaseState(tracker, compTrack1);
            CallableLoci.CallableBaseState comp2 = getCallableBaseState(tracker, compTrack2);

            if ( printState != null && comp1.getState() == printState1 && comp2.getState() == printState2 ) {
                out.printf("%s %s %s %s%n", comp1.getLocation(), comp1.getState(), comp2.getLocation(), comp2.getState());
            }

            return Arrays.asList(comp1, comp2);
        } else {
            return null;
        }
    }

    private CallableLoci.CallableBaseState getCallableBaseState(RefMetaDataTracker tracker, RodBinding<BEDFeature> rodBinding) {
        //System.out.printf("tracker %s%n", tracker);
        List<BEDFeature> bindings = tracker.getValues(rodBinding);
        if ( bindings.size() != 1 ) {
            throw new UserException.MalformedFile(String.format("%s track isn't a properly formatted CallableBases object!", rodBinding.getName()));
        }

        BEDFeature bed = bindings.get(0);
        GenomeLoc loc = getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), bed.getStart(), bed.getEnd());
        CallableLoci.CalledState state = CallableLoci.CalledState.valueOf(bed.getName());
        return new CallableLoci.CallableBaseState(getToolkit().getGenomeLocParser(),loc, state);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    // --------------------------------------------------------------------------------------------------------------
    public long[][] reduceInit() {
        int n = CallableLoci.CalledState.values().length;
        return new long[n][n];
    }

    public long[][] reduce(List<CallableLoci.CallableBaseState> comps, long[][] sum) {
        if ( comps != null ) {
            CallableLoci.CallableBaseState comp1 = comps.get(0);
            CallableLoci.CallableBaseState comp2 = comps.get(1);

            sum[comp1.getState().ordinal()][comp2.getState().ordinal()]++;
        }

        return sum;
    }

    public void onTraversalDone(long[][] result) {
        for ( CallableLoci.CalledState state1 : CallableLoci.CalledState.values() ) {
            for ( CallableLoci.CalledState state2 : CallableLoci.CalledState.values() ) {
                out.printf("%s %s %s %s %d%n", compTrack1.getName(), compTrack2.getName(), state1, state2, result[state1.ordinal()][state2.ordinal()]);
            }
        }
    }
}