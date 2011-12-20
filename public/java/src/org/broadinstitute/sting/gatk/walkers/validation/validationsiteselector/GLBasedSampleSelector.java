/*
 * Copyright (c) 2010, The Broad Institute
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
package org.broadinstitute.sting.gatk.walkers.validation.validationsiteselector;

import org.broadinstitute.sting.gatk.walkers.genotyper.AlleleFrequencyCalculationResult;
import org.broadinstitute.sting.gatk.walkers.genotyper.ExactAFCalculationModel;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;


public class GLBasedSampleSelector extends SampleSelector {
    Map<Integer,double[][]> numAllelePriorMatrix = new HashMap<Integer,double[][]>();
    double referenceLikelihood;
    public GLBasedSampleSelector(TreeSet<String> sm, double refLik) {
        super(sm);
        referenceLikelihood = refLik;
    }

    public boolean selectSiteInSamples(VariantContext vc) {
        if ( samples == null || samples.isEmpty() )
            return true;
        // want to include a site in the given samples if it is *likely* to be variant (via the EXACT model)
        // first subset to the samples
        VariantContext subContext = vc.subContextFromSamples(samples);

        // now check to see (using EXACT model) whether this should be variant
        // do we want to apply a prior? maybe user-spec?
        double[][] flatPrior = createFlatPrior(vc.getAlleles());
        AlleleFrequencyCalculationResult result = new AlleleFrequencyCalculationResult(vc.getAlternateAlleles().size(),2*samples.size());
        ExactAFCalculationModel.linearExactMultiAllelic(subContext.getGenotypes(),vc.getAlternateAlleles().size(),flatPrior,result,true);
        // do we want to let this qual go up or down?
        if ( result.getLog10PosteriorOfAFzero() < referenceLikelihood ) {
            return true;
        }

        return false;
    }

    private double[][] createFlatPrior(List<Allele> alleles) {
        if ( ! numAllelePriorMatrix.containsKey(alleles.size()) ) {
            numAllelePriorMatrix.put(alleles.size(), new double[alleles.size()][1+2*samples.size()]);
        }

        return numAllelePriorMatrix.get(alleles.size());
    }
}
