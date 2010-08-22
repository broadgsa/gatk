/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.EnumSet;
import java.io.File;
import java.io.PrintStream;

/**
 * Test routine for new VariantContext object
 */
public class TestVariantContextWalker extends RodWalker<Integer, Integer> {
    @Output
    PrintStream out;

    @Argument(fullName="takeFirstOnly", doc="Only take the first second at a locus, as opposed to all", required=false)
    boolean takeFirstOnly = false;

    @Argument(fullName="onlyContextsOfType", doc="Only take variant contexts of this type", required=false)
    VariantContext.Type onlyOfThisType = null;

    @Argument(fullName="onlyContextsStartinAtCurrentPosition", doc="Only take variant contexts at actually start at the current position, excluding those at span to the current location but start earlier", required=false)
    boolean onlyContextsStartinAtCurrentPosition = false;

    @Argument(fullName="printPerLocus", doc="If true, we'll psetenv LD_LIBRARY_PATH .:/util/gcc-4.3.0/lib64:/util/gcc-4.3.0/lib/gcc/x86_64-unknown-linux-gnu/4.3.0:/util/gcc-4.3.0/lib/:/util/lib:/broad/tools/Linux/x86_64/pkgs/boost_1.38.0/lib:/humgen/gsa-scr1/GATK_Data/bwarint the variant contexts, in addition to counts", required=false)
    boolean printContexts = false;

    @Argument(fullName="outputVCF", doc="If provided, we'll convert the first input context into a VCF", required=false)
    String outputVCF = null;

    private VCFWriter writer = null;
    private boolean wroteHeader = false;

    public void initialize() {
        if ( outputVCF != null )
            writer = new StandardVCFWriter(new File(outputVCF));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( ref == null )
            return 0;
        else {
            EnumSet<VariantContext.Type> allowedTypes = onlyOfThisType == null ? null : EnumSet.of(onlyOfThisType);

            int n = 0;
            for (VariantContext vc : tracker.getAllVariantContexts(ref, allowedTypes, context.getLocation(), onlyContextsStartinAtCurrentPosition, takeFirstOnly) ) {
                if ( writer != null && n == 0 ) {
                    if ( ! wroteHeader ) {
                        writer.writeHeader(VariantContextAdaptors.createVCFHeader(null, vc));
                        wroteHeader = true;
                    }

                    writer.add(vc, ref.getBase());
                }

                n++;
                if ( printContexts ) out.printf("       %s%n", vc);
            }

            if ( n > 0 && printContexts ) {
                out.printf("%s => had %d variant context objects%n", context.getLocation(), n);
                out.printf("---------------------------------------------%n");
            }
            
            return n;
        }
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer point, Integer sum) {
        return point + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
        // Double check traversal result to make count is the same.
        // TODO: Is this check necessary?
        out.println("[REDUCE RESULT] Traversal result is: " + result);
    }
}