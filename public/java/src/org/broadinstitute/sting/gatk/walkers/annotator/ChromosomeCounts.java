/*
 * Copyright (c) 2010.
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

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;


/**
 * Allele count in genotypes, for each ALT allele, in the same order as listed;
 * allele Frequency, for each ALT allele, in the same order as listed; total number
 * of alleles in called genotypes.
 */
public class ChromosomeCounts extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    public static final String[] keyNames = { VCFConstants.ALLELE_NUMBER_KEY, VCFConstants.ALLELE_COUNT_KEY, VCFConstants.ALLELE_FREQUENCY_KEY };
    public static final VCFInfoHeaderLine[] descriptions = {
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY),
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY),
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY) };

    private Set<String> founderIds = new HashSet<String>();

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( ! vc.hasGenotypes() )
            return null;

        return VariantContextUtils.calculateChromosomeCounts(vc, new HashMap<String, Object>(), true,founderIds);
    }

    public void initialize ( AnnotatorCompatibleWalker walker, GenomeAnalysisEngine toolkit, Set<VCFHeaderLine> headerLines ){
        //If families were given, get the founders ids
        founderIds = ((Walker)walker).getSampleDB().getFounderIds();
    }

    public Map<String, Object> annotate(Map<String, Map<Allele, List<GATKSAMRecord>>> stratifiedContexts, VariantContext vc) {
        if ( ! vc.hasGenotypes() )
            return null;

        return VariantContextUtils.calculateChromosomeCounts(vc, new HashMap<String, Object>(), true);
    }

    public List<String> getKeyNames() {
        return Arrays.asList(keyNames);
    }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(descriptions); }
}