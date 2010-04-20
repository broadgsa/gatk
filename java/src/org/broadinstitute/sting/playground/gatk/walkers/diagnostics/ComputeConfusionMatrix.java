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

package org.broadinstitute.sting.playground.gatk.walkers.diagnostics;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.BaseUtils;
import net.sf.samtools.SAMRecord;

import java.util.HashMap;
import java.util.Arrays;
import java.util.Hashtable;

/**
 * Computes empirical base confusion matrix, and optionally computes
 * these matrices with up to five bases of preceding context
 */
@Reference(window=@Window(start=-5,stop=5))
public class ComputeConfusionMatrix extends LocusWalker<Integer, Integer> {
    @Argument(fullName="minimumDepth", shortName="minDepth", doc="Require locus pileup to have specified minimum depth", required=false)
    public Integer MIN_DEPTH = 10;

    @Argument(fullName="maximumDepth", shortName="maxDepth", doc="Require locus pileup to have specified maximum depth", required=false)
    public Integer MAX_DEPTH = 100;

    @Argument(fullName="contextWindowSize", shortName="window", doc="Size of context window", required=false)
    public Integer WINDOW_SIZE = 0;

    private Hashtable<String, Integer> confusionCounts = new Hashtable<String, Integer>();

    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int pileupSize = context.size();

        int numAlts = 0;
        int[] baseCounts = context.getBasePileup().getBaseCounts();
        for (int baseIndex = 0; baseIndex < baseCounts.length; baseIndex++) {
            if (baseIndex != ref.getBaseIndex()) {
                numAlts += baseCounts[baseIndex];
            }
        }

        return (
            pileupSize >= MIN_DEPTH &&        // don't process regions without a reasonable pileup
            pileupSize < MAX_DEPTH &&         // don't process suspiciously overcovered regions
            ref.getBases().length % 2 == 1 && // don't process regions that don't have a full context window
            numAlts == 1 &&                   // don't process regions that have more than one mismatching base
            ref.getBaseIndex() >= 0           // don't process a locus with an ambiguous reference base
        );
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int windowLength = ref.getBases().length;
        int windowCenter = (windowLength - 1)/2;

        String fwRefBases = new String(ref.getBases());
        String fwRefBase = String.format("%c", ref.getBase());
        String fwWindowLeft = fwRefBases.substring(windowCenter - WINDOW_SIZE, windowCenter);

        //String rcRefBases = new String(BaseUtils.simpleReverseComplement(ref.getBases()));
        //String rcRefBase = String.format("%c", BaseUtils.simpleComplement(ref.getBase()));
        //String rcWindowRight = rcRefBases.substring(windowCenter + 1, windowCenter + 1 + WINDOW_SIZE);

        int[] baseCounts = context.getBasePileup().getBaseCounts();
        int altBaseIndex = -1;
        for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
            if (baseCounts[baseIndex] == 1) {
                altBaseIndex = baseIndex;
            }
        }
        
        String fwAltBase = String.format("%c", BaseUtils.baseIndexToSimpleBase(altBaseIndex));
        //String rcAltBase = BaseUtils.simpleComplement(fwAltBase);

        for (int readIndex = 0; readIndex < context.getReads().size(); readIndex++) {
            SAMRecord read = context.getReads().get(readIndex);
            int offset = context.getOffsets().get(readIndex);

            char base = read.getReadString().charAt(offset);
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);

            if (baseIndex == altBaseIndex) {
                if (read.getReadNegativeStrandFlag()) {
                    //incrementConfusionCounts(rcWindowRight, rcRefBase, rcAltBase);
                } else {
                    incrementConfusionCounts(fwWindowLeft, fwAltBase, fwRefBase);
                }
            }
        }

        return null;
    }

    private void incrementConfusionCounts(String context, String altBase, String refBase) {
        String key = String.format("%s:%s:%s", context, altBase, refBase);

        Integer counts = confusionCounts.get(key);
        if (counts == null) { counts = 0; }

        confusionCounts.put(key, counts + 1);
    }

    public Integer reduceInit() {
        return null;
    }

    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer result) {
        String[] keys = confusionCounts.keySet().toArray(new String[0]);
        Arrays.sort(keys);

        HashMap<String, Integer> contextualNorms = new HashMap<String, Integer>();
        for (String key : keys) {
            String[] fields = key.split(":");

            String contextualKey = String.format("%s:%s", fields[0], fields[1]);
            Integer contextualCount = contextualNorms.get(contextualKey);
            if (contextualCount == null) { contextualCount = 0; }
            contextualNorms.put(contextualKey, contextualCount + confusionCounts.get(key));
        }

        out.printf("confusionMatrix\tcontext\talt\tref\tcounts\ttotal\tfraction\n");
        for (String key : keys) {
            String[] fields = key.split(":");
            String contextualKey = String.format("%s:%s", fields[0], fields[1]);

            out.printf(
                "confusionMatrix\t%s\t%s\t%s\t%d\t%d\t%f\n",
                fields[0],
                fields[1],
                fields[2],
                confusionCounts.get(key),
                contextualNorms.get(contextualKey),
                confusionCounts.get(key)/((float) contextualNorms.get(contextualKey))
            );
        }
    }
}
