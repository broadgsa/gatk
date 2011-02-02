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

package org.broadinstitute.sting.oneoffprojects.walkers.phasing;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.phasing.ReadBackedPhasingWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Emits specific fields as dictated by the user from one or more VCF files.
 */
@Requires(value = {})
public class PhasingEval extends RodWalker<Integer, Integer> {

    @Output(doc = "File to which results should be written", required = true)
    protected PrintStream out;

    @Argument(doc = "sample to emit", required = false)
    protected String sample = null;

    @Argument(doc = "Analysis to perform", required = true)
    protected Analysis analysis;

    public enum Analysis {
        PHASING_BY_AC
    }

    private class PhasingByAC {
        int myAC = 0;
        int myAN = 0;
        int nHets = 0;
        int nHetsPhased = 0;

        public PhasingByAC(int myAC, int myAN) {
            this.myAC = myAC;
            this.myAN = myAN;
        }
    }

    List<PhasingByAC> phasingByACs = new ArrayList<PhasingByAC>();

    public void initialize() {
        Set<String> samples = SampleUtils.getSampleList(VCFUtils.getVCFHeadersFromRods(getToolkit(), null));
        int AN = 2 * samples.size();
        for (int i = 0; i <= AN; i++) {
            phasingByACs.add(new PhasingByAC(i, AN));
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null) // RodWalkers can make funky map calls
            return 0;

        Collection<VariantContext> vcs = tracker.getAllVariantContexts(ref, context.getLocation());
        for (VariantContext vc : vcs) {
            if (sample != null)
                vc = vc.subContextFromGenotypes(vc.getGenotype(sample));

            if (analysis == Analysis.PHASING_BY_AC) {
                int homref = vc.getHomRefCount();
                int homalt = vc.getHomVarCount();
                int het = vc.getHetCount();
                int ac = 2 * homalt + het;

                //int an = 2 * (homref + homalt + het);

                PhasingByAC data = phasingByACs.get(ac);
                data.nHets += het > 0 ? 1 : 0;
                data.nHetsPhased += isPhysicallyPhased(vc.getGenotypes().values()) ? 1 : 0;
            }
        }

        return 1;
    }

    private boolean isPhysicallyPhased(Collection<Genotype> genotypes) {
        for (Genotype g : genotypes) {
            if (g.isHet() && g.hasAttribute(ReadBackedPhasingWalker.PQ_KEY))
                return true;
        }

        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {
        if (analysis == Analysis.PHASING_BY_AC) {
            out.println(Utils.join("\t", Arrays.asList("ac", "an", "nhets", "nhetphased")));
            for (PhasingByAC pac : phasingByACs) {
                out.printf("%d\t%d\t%d\t%d%n", pac.myAC, pac.myAN, pac.nHets, pac.nHetsPhased);
            }
        }
    }
}