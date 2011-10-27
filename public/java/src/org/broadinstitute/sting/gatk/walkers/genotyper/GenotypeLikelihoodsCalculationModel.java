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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Map;


/**
 * The model representing how we calculate genotype likelihoods
 */
public abstract class GenotypeLikelihoodsCalculationModel implements Cloneable {

    public enum Model {
        SNP,
        INDEL,
        BOTH
    }

    public enum GENOTYPING_MODE {
        /** the Unified Genotyper will choose the most likely alternate allele */
        DISCOVERY,
        /** only the alleles passed in from a VCF rod bound to the -alleles argument will be used for genotyping */
        GENOTYPE_GIVEN_ALLELES
    }

    protected UnifiedArgumentCollection UAC;
    protected Logger logger;

    /**
     * Create a new object
     * @param logger        logger
     * @param UAC           unified arg collection
     */
    protected GenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        if ( logger == null || UAC == null ) throw new ReviewedStingException("Bad arguments");
        this.UAC = UAC.clone();
        this.logger = logger;
    }

    /**
     * Must be overridden by concrete subclasses
     *
     * @param tracker              rod data
     * @param ref                  reference context
     * @param contexts             stratified alignment contexts
     * @param contextType          stratified context type
     * @param priors               priors to use for GLs
     * @param GLs                  hash of sample->GL to fill in
     * @param alternateAlleleToUse the alternate allele to use, null if not set
     *
     * @param useBAQedPileup
     * @return genotype likelihoods per sample for AA, AB, BB
     */
    public abstract Allele getLikelihoods(RefMetaDataTracker tracker,
                                          ReferenceContext ref,
                                          Map<String, AlignmentContext> contexts,
                                          AlignmentContextUtils.ReadOrientation contextType,
                                          GenotypePriors priors,
                                          Map<String, MultiallelicGenotypeLikelihoods> GLs,
                                          Allele alternateAlleleToUse, boolean useBAQedPileup);

    protected int getFilteredDepth(ReadBackedPileup pileup) {
        int count = 0;
        for ( PileupElement p : pileup ) {
            if ( BaseUtils.isRegularBase( p.getBase() ) )
                count++;
        }

        return count;
    }
}