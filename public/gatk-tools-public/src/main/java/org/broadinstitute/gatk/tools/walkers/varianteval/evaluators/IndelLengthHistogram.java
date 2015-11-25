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

package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Molten;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Simple utility for histogramming indel lengths
 *
 * Based on code from chartl
 *
 * @author Mark DePristo
 * @since 3/21/12
 */
@Analysis(description = "Indel length histogram", molten = true)
public class IndelLengthHistogram extends VariantEvaluator implements StandardEval {
    private final Map<Integer, Integer> counts = new HashMap<Integer, Integer>();
    private final static boolean asFrequencies = true;
    int nIndels = 0;

    @Molten(variableName = "Length", valueName = "Freq", variableFormat = "%d", valueFormat = "%.2f")
    public TreeMap<Object, Object> results;
    
    public final static int MAX_SIZE_FOR_HISTOGRAM = 10;
    private final static boolean INCLUDE_LONG_EVENTS_AT_MAX_SIZE = false;

    public IndelLengthHistogram() {
        initializeCounts(MAX_SIZE_FOR_HISTOGRAM);
    }

    private void initializeCounts(int size) {
        for ( int i = -size; i <= size; i++ ) {
            if ( i != 0 ) counts.put(i, 0);
        }
    }

    @Override
    public void finalizeEvaluation() {
        if ( asFrequencies ) {
            results = new TreeMap<Object, Object>();
            for ( final int len : counts.keySet() ) {
                final double value = nIndels == 0 ? 0.0 : counts.get(len) / (1.0 * nIndels);
                results.put(len, value);
            }
        } else {
            results = new TreeMap<Object, Object>(results);
        }
    }

    @Override
    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void update1(final VariantContext eval, final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( eval.isIndel() && ! eval.isComplexIndel() ) {
            if ( ! ( getWalker().ignoreAC0Sites() && eval.isMonomorphicInSamples() )) {
                // only if we are actually polymorphic in the subsetted samples should we count the allele
                for ( Allele alt : eval.getAlternateAlleles() ) {
                    final int alleleSize = alt.length() - eval.getReference().length();
                    if ( alleleSize == 0 ) throw new ReviewedGATKException("Allele size not expected to be zero for indel: alt = " + alt + " ref = " + eval.getReference());
                    updateLengthHistogram(eval.getReference(), alt);
                }
            }
        }
    }

    /**
     * Update the histogram with the implied length of the indel allele between ref and alt (alt.len - ref.len).
     *
     * If this size is outside of MAX_SIZE_FOR_HISTOGRAM, the size is capped to MAX_SIZE_FOR_HISTOGRAM,
     * if INCLUDE_LONG_EVENTS_AT_MAX_SIZE is set.
     *
     * @param ref
     * @param alt
     */
    public void updateLengthHistogram(final Allele ref, final Allele alt) {
        int len = alt.length() - ref.length();
        if ( INCLUDE_LONG_EVENTS_AT_MAX_SIZE ) {
            if ( len > MAX_SIZE_FOR_HISTOGRAM ) len = MAX_SIZE_FOR_HISTOGRAM;
            if ( len < -MAX_SIZE_FOR_HISTOGRAM ) len = -MAX_SIZE_FOR_HISTOGRAM;
        }
        
        if ( Math.abs(len) > MAX_SIZE_FOR_HISTOGRAM )
            return;
        
        nIndels++;
        counts.put(len, counts.get(len) + 1);
    }
}
