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

package org.broadinstitute.gatk.tools.walkers.annotator.interfaces;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.tools.walkers.annotator.ReducibleAnnotationData;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.List;
import java.util.Map;

/**
 * An interface for annotations that are calculated using raw data across samples, rather than the median (or median of median) of samples values
 */
public interface ReducibleAnnotation extends AnnotationType {
    public abstract String getRawKeyName();

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     *
     * @param tracker
     * @param walker
     * @param ref
     * @param stratifiedContexts
     * @param vc
     * @param stratifiedPerReadAlleleLikelihoodMap
     * @return
     */
    public abstract Map<String, Object> annotateRawData(final RefMetaDataTracker tracker,
                                                        final AnnotatorCompatible walker,
                                                        final ReferenceContext ref,
                                                        final Map<String, AlignmentContext> stratifiedContexts,
                                                        final VariantContext vc,
                                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap);

    /**
     * Combine raw data, typically during the merging of raw data contained in multiple gVCFs as in CombineGVCFs and the
     * preliminary merge for GenotypeGVCFs
     * @param allelesList   The merged allele list across all variants being combined/merged
     * @param listOfRawData The raw data for all the variants being combined/merged
     * @return
     */
    public abstract Map<String, Object> combineRawData(final List<Allele> allelesList, final List <? extends ReducibleAnnotationData> listOfRawData);


    /**
     * Calculate the final annotation value from the raw data
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    public abstract Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC);

    /**
     *
     * @param vc
     * @param pralm
     * @param rawAnnotations
     */
    public abstract void calculateRawData(VariantContext vc, Map<String, PerReadAlleleLikelihoodMap> pralm, ReducibleAnnotationData rawAnnotations);
}
