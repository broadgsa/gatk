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

package org.broadinstitute.gatk.tools.walkers.annotator;

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.BaseUtils;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Percentage of N bases
 *
 * <p>N occurs in a sequence when the sequencer does not have enough information to determine which base it should call. The presence of many Ns at the same site lowers our confidence in any calls made there, because it suggests that there was some kind of technical difficulty that interfered with the sequencing process.</p>
 *
 * <h3>Note</h3>
 * <p>In GATK versions 3.2 and earlier, this annotation only counted N bases from reads generated with SOLiD technology. This functionality was generalized for all sequencing platforms in GATK version 3.3.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_BaseCounts.php">BaseCounts</a></b> counts the number of A, C, G, T bases across all samples.</li>
 * </ul>
 *
 * */
public class NBaseCount extends InfoFieldAnnotation {
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if( stratifiedContexts.size() == 0 )
            return null;

        int countNBase = 0;
        int countRegularBase = 0;

        for( final AlignmentContext context : stratifiedContexts.values() ) {
            for( final PileupElement p : context.getBasePileup()) {
                if( BaseUtils.isNBase( p.getBase() ) ) {
                    countNBase++;
                } else if( BaseUtils.isRegularBase( p.getBase() ) ) {
                    countRegularBase++;
                }
            }
        }
        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.format("%.4f", (double)countNBase / (double)(countNBase + countRegularBase + 1)));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.N_BASE_COUNT_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0))); }
}
