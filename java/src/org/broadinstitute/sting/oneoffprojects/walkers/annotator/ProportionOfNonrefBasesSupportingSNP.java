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
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;

import java.util.Map;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 17, 2009
 * Time: 2:48:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProportionOfNonrefBasesSupportingSNP implements InfoFieldAnnotation {
    private String KEY_NAME = "prop_nonref_that_are_snp";

    public String getKeyName() { return KEY_NAME; }

    public VCFInfoHeaderLine getDescription() {
        return new VCFInfoHeaderLine(KEY_NAME,
                        1,VCFInfoHeaderLine.INFO_TYPE.Float,"Simple proportion of non-reference bases that are the SNP base");
    }

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> context, VariantContext vc) {
        if ( ! vc.isSNP() || ! vc.isBiallelic() )
            return null;

        Pair<Integer,Integer> totalNonref_totalSNP = new Pair<Integer,Integer>(0,0);
        for ( String sample : context.keySet() ) {
            ReadBackedPileup pileup = context.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            totalNonref_totalSNP = getNonrefAndSNP(pileup, ref.getBaseAsChar(), vc.getAlternateAllele(0).toString().charAt(0), totalNonref_totalSNP);

        }
        if ( totalNonref_totalSNP.equals(new Pair<Integer,Integer>(0,0)) )
            return null;
        double p = getProportionOfNonrefBasesThatAreSNP(totalNonref_totalSNP);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%f", p ));
        return map;
    }

    private Pair<Integer,Integer> getNonrefAndSNP(ReadBackedPileup p, char ref, char snp, Pair<Integer,Integer> totals) {
        int[] counts = p.getBaseCounts();
        int nonrefCounts = 0;
        int snpCounts = counts[BaseUtils.simpleBaseToBaseIndex(snp)];
        for ( char c : BaseUtils.BASES ) {
            if ( ! BaseUtils.basesAreEqual((byte) c, (byte) ref) ) {
                nonrefCounts += counts[BaseUtils.simpleBaseToBaseIndex(c)];
            }
        }

        totals.first+=nonrefCounts;
        totals.second+=snpCounts;
        return totals;
    }

    private double getProportionOfNonrefBasesThatAreSNP( Pair<Integer,Integer> totalNonref_totalSNP ) {
        return ( 1.0 + totalNonref_totalSNP.second ) / (1.0 + totalNonref_totalSNP.first );
    }
}
