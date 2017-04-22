/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.Collectors;


/**
 * Count of A, C, G, T bases across all samples
 *
 * <p> This annotation returns the counts of A, C, G, and T bases across all samples, in that order.</p>
 * <h3>Example:</h3>
 *
 * <pre>BaseCounts=3,0,3,0</pre>
 *
 * <p>
 *     This means the number of A bases seen is 3, the number of T bases seen is 0, the number of G bases seen is 3, and the number of T bases seen is 0.
 * </p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_NBaseCount.php">NBaseCount</a></b> counts the percentage of N bases.</li>
 * </ul>
 */

 public class BaseCounts extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {

        if ( stratifiedPerReadAlleleLikelihoodMap == null || stratifiedPerReadAlleleLikelihoodMap.isEmpty() ) {
            return null;
        }

        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), Arrays.stream(getBaseCounts(stratifiedPerReadAlleleLikelihoodMap, vc)).boxed().collect(Collectors.toList()));
        return map;
    }

    /**
     * Counts of observed bases at a genomic position (e.g. {13,0,0,1} at chr1:100,000,000) over all samples
     *
     * @param stratifiedPerReadAlleleLikelihoodMap for each read, the underlying alleles represented by an aligned read, and corresponding relative likelihood.
     * @param vc variant context
     * @return count of A, C, G, T bases
     */
    private int[] getBaseCounts(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap, final VariantContext vc) {
        final Set<Allele> alleles = new HashSet<>(vc.getAlleles());
        final int[] baseCounts = new int[4];
        for ( final Map.Entry<String, PerReadAlleleLikelihoodMap> strat : stratifiedPerReadAlleleLikelihoodMap.entrySet() ) {
            final int[] counts = AlignmentUtils.countBasesAtPileupPosition(strat.getValue(), alleles, vc.getStart());
            for ( int i = 0; i < baseCounts.length; i++ ) {
                baseCounts[i] += counts[i];
            }
        }

        return baseCounts;
    }

    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.BASE_COUNTS_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0))); }
}