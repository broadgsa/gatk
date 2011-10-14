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

package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

/**
 * Test routine for new VariantContext object
 */
@Reference(window=@Window(start=-20,stop=1))
public class TestVariantContextWalker extends RodWalker<Integer, Integer> {
    @Output
    PrintStream out;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(fullName="takeFirstOnly", doc="Only take the first second at a locus, as opposed to all", required=false)
    boolean takeFirstOnly = false;

    @Argument(fullName="onlyContextsOfType", doc="Only take variant contexts of this type", required=false)
    VariantContext.Type onlyOfThisType = null;

    @Argument(fullName="onlyContextsStartinAtCurrentPosition", doc="Only take variant contexts at actually start at the current position, excluding those at span to the current location but start earlier", required=false)
    boolean onlyContextsStartinAtCurrentPosition = false;

    @Argument(fullName="printPerLocus", doc="If true, we'll print the variant contexts, in addition to counts", required=false)
    boolean printContexts = false;

    @Argument(fullName="outputVCF", doc="If provided, we'll convert the first input context into a VCF", required=false)
    VCFWriter writer = null;

    private boolean wroteHeader = false;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( ref == null )
            return 0;
        else {
            EnumSet<VariantContext.Type> allowedTypes = onlyOfThisType == null ? null : EnumSet.of(onlyOfThisType);

            int n = 0;
            List<VariantContext> contexts;
            if ( onlyContextsStartinAtCurrentPosition )
                contexts = tracker.getValues(variantCollection.variants, context.getLocation());
            else // ! onlyContextsStartinAtCurrentPosition
                contexts = tracker.getValues(variantCollection.variants);

            for ( VariantContext vc : contexts ) {
                if ( allowedTypes == null || allowedTypes.contains(vc.getType()) ) {
                    // we need to trigger decoding of the genotype string to pass integration tests
                    vc.getGenotypes();

                    if ( writer != null && n == 0 ) {
                        if ( ! wroteHeader ) {
                            writer.writeHeader(VariantContextAdaptors.createVCFHeader(null, vc));
                            wroteHeader = true;
                        }

                        writer.add(vc);
                    }

                    n++;
                    if ( printContexts ) out.printf("       %s%n", vc);
                    if ( takeFirstOnly ) break;
                }
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
