/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class TandemRepeatAnnotator extends InfoFieldAnnotation implements StandardAnnotation {
    private static final String STR_PRESENT = "STR";
    private static final String REPEAT_UNIT_KEY = "RU";
    private static final String REPEATS_PER_ALLELE_KEY = "RPA";
    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( !vc.isIndel())
            return null;

        Pair<List<Integer>,byte[]> result = VariantContextUtils.getNumTandemRepeatUnits(vc, ref.getForwardBases());
        if (result == null)
            return null;

        byte[] repeatUnit = result.second;
        List<Integer> numUnits = result.first;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(STR_PRESENT,true);
        map.put(REPEAT_UNIT_KEY,new String(repeatUnit));
        map.put(REPEATS_PER_ALLELE_KEY, numUnits);

        return map;
    }

    public static final String[] keyNames = {STR_PRESENT, REPEAT_UNIT_KEY,REPEATS_PER_ALLELE_KEY };
    public static final VCFInfoHeaderLine[] descriptions = {
            new VCFInfoHeaderLine(STR_PRESENT, 0, VCFHeaderLineType.Flag, "Variant is a short tandem repeat"),
            new VCFInfoHeaderLine(REPEAT_UNIT_KEY, 1, VCFHeaderLineType.String, "Tandem repeat unit (bases)"),
            new VCFInfoHeaderLine(REPEATS_PER_ALLELE_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Number of times tandem repeat unit is repeated, for each allele (including reference)") };

    public List<String> getKeyNames() {
        return Arrays.asList(keyNames);
    }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(descriptions); }

}
