/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.variantutils.VariantsToTable;
import org.broadinstitute.sting.utils.Utils;

import java.io.PrintStream;
import java.util.*;

/**
 * Emits specific fields as dictated by the user from one or more VCF files.
 */
@Requires(value={})
public class AssessSimulatedPerformance extends RodWalker<Integer, Integer> {
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(fullName="fields", shortName="F", doc="Fields to emit from the calls VCF", required=false)
    public String FIELDS = "CHROM,POS,REF,ALT,QUAL,AC,AN,DP,Q,MODE";

    @Argument(fullName="maxRecords", shortName="M", doc="Maximum number of records to emit, if provided", required=false)
    public int MAX_RECORDS = -1;
    int nRecords = 0;

    private List<String> fieldsToTake;

    public void initialize() {
        fieldsToTake = Arrays.asList(FIELDS.split(","));

        for ( String source : Arrays.asList("sim", "called")) {
            out.print(source + "." + Utils.join("\t" + source + ".", fieldsToTake));
            out.print("\t");
        }
        out.println();
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        if ( ++nRecords < MAX_RECORDS || MAX_RECORDS == -1 ) {
            printVCFields("sim", tracker, ref, context);
            printVCFields("called", tracker, ref, context);
            out.println();
            return 1;
        } else {
            if ( nRecords >= MAX_RECORDS ) {
                logger.warn("Calling sys exit to leave after " + nRecords + " records");
                System.exit(0); // todo -- what's the recommend way to abort like this?
            }
            return 0;
        }
    }

    private void printVCFields(String name, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariantContext vc = tracker.getVariantContext(ref, name, null, context.getLocation(), true);
        out.print(Utils.join("\t", VariantsToTable.extractFields(vc, fieldsToTake)));
        out.print("\t");
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {}
}