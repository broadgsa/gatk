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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;


public class AlleleBalance implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        if ( !vc.isBiallelic() )
            return null;
        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( !vc.hasGenotypes() )
            return null;

        double ratio = 0.0;
        double totalWeights = 0.0;
        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {
            // we care only about het calls
            if ( !genotype.getValue().isHet() )
                continue;

            StratifiedAlignmentContext context = stratifiedContexts.get(genotype.getKey());
            if ( context == null )
                continue;

            if ( vc.isSNP() ) {
                final String bases = new String(context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup().getBases());
                if ( bases.length() == 0 )
                    return null;
                char refChr = vc.getReference().toString().charAt(0);
                char altChr = vc.getAlternateAllele(0).toString().charAt(0);

                int refCount = MathUtils.countOccurrences(refChr, bases);
                int altCount = MathUtils.countOccurrences(altChr, bases);

                // sanity check
                if ( refCount + altCount == 0 )
                    continue;

                // weight the allele balance by genotype quality so that e.g. mis-called homs don't affect the ratio too much
                ratio += genotype.getValue().getNegLog10PError() * ((double)refCount / (double)(refCount + altCount));
                totalWeights += genotype.getValue().getNegLog10PError();
            } else if ( vc.isIndel() ) {
                final ReadBackedExtendedEventPileup indelPileup = context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getExtendedEventPileup();
                if ( indelPileup == null ) {
                    continue;
                }
                // todo -- actually care about indel length from the pileup (agnostic at the moment)
                int refCount = indelPileup.size();
                int altCount = vc.isInsertion() ? indelPileup.getNumberOfInsertions() : indelPileup.getNumberOfDeletions();

                if ( refCount + altCount == 0 ) {
                    continue;
                }

                ratio += /* todo -- make not uniform */ 1 * ((double) refCount) / (double) (refCount + altCount);
                totalWeights += 1;
            }
        }

        // make sure we had a het genotype
        if ( MathUtils.compareDoubles(totalWeights, 0.0) == 0 )
            return null;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", (ratio / totalWeights)));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("AB"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("AB", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Allele Balance for hets (ref/(ref+alt))")); }
}
