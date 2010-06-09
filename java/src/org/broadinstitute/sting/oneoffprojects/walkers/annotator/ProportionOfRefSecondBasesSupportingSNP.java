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

package org.broadinstitute.sting.oneoffprojects.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 17, 2009
 * Time: 2:18:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProportionOfRefSecondBasesSupportingSNP implements InfoFieldAnnotation {
    private String KEY_NAME = "ref_2bb_snp_prop";
    private boolean USE_MAPQ0_READS = false;

    public List<String> getKeyNames() { return Arrays.asList(KEY_NAME); }

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> context, VariantContext vc) {
        if ( ! vc.isSNP() || ! vc.isBiallelic() )
            return null;

        Pair<Integer,Integer> totalAndSNPSupporting = new Pair<Integer,Integer>(0,0);
        for ( String sample : context.keySet() ) {
            ReadBackedPileup pileup = context.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            totalAndSNPSupporting = getTotalRefAndSNPSupportCounts(pileup, ref.getBaseAsChar(), vc.getAlternateAllele(0).toString().charAt(0), totalAndSNPSupporting);

        }
        if ( totalAndSNPSupporting.equals(new Pair<Integer,Integer>(0,0)) )
            return null;
        
        double p = getProportionOfRefSecondaryBasesSupportingSNP(totalAndSNPSupporting);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%f", p ));
        return map;
    }

    private double getProportionOfRefSecondaryBasesSupportingSNP(Pair<Integer,Integer> tRef_snpSupport) {
        return ( 1.0 + tRef_snpSupport.second) / (1.0 + tRef_snpSupport.first );
    }

    private Pair<Integer,Integer> getTotalRefAndSNPSupportCounts(ReadBackedPileup p, char ref, char snp, Pair<Integer,Integer> refAndSNPCounts) {
        int nRefBases = 0;
        int nSecondBasesSupportingSNP = 0;
        for (PileupElement e : p ) {
            if ( BaseUtils.basesAreEqual( e.getBase(), (byte) ref ) ) {
                if ( BaseUtils.isRegularBase(e.getSecondBase()) ) {
                    nRefBases++;
                    if ( BaseUtils.basesAreEqual( e.getSecondBase(), (byte) snp ) ) {
                        nSecondBasesSupportingSNP++;
                    }
                }
            }
        }

        refAndSNPCounts.first+=nRefBases;
        refAndSNPCounts.second+=nSecondBasesSupportingSNP;
        return refAndSNPCounts;
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(KEY_NAME,
                        1,VCFInfoHeaderLine.INFO_TYPE.Float,"Simple proportion of second best base calls for reference base that support the SNP base"));
    }
}
