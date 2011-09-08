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

package org.broadinstitute.sting.gatk.walkers.cancer;

import net.sf.picard.util.MathUtil;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Assigns somatic status to a set of calls
 */
public class AssignSomaticStatus extends RodWalker<Integer, Integer> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(shortName="t", fullName="tumorsample", required=true, doc="List of tumor samples")
    public Set<String> tumorSamplesArg;

    @Argument(shortName="somaticPriorQ", fullName="somaticPriorQ", required=false, doc="Phred-scaled probability that a site is a somatic mutation")
    public byte somaticPriorQ = 60;

    @Output
    protected VCFWriter vcfWriter = null;

    private final String SOMATIC_TAG_NAME = "SOMATIC";
    private final String SOURCE_NAME = "AssignSomaticStatus";

    private Set<String> tumorSamples = new HashSet<String>();
    private Set<String> normalSamples = new HashSet<String>();

    /**
     * Parse the familial relationship specification, and initialize VCF writer
     */
    public void initialize() {
        List<String> rodNames = new ArrayList<String>();
        rodNames.add(variantCollection.variants.getName());

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        // set up tumor and normal samples
        for ( final String sample : vcfSamples ) {
            if ( tumorSamplesArg.contains(sample) )
                tumorSamples.add(sample);
            else
                normalSamples.add(sample);
        }
        logger.info("N tumor  samples: " + tumorSamples.size());
        logger.info("N normal samples: " + normalSamples.size());
        if ( tumorSamples.size() != normalSamples.size() )
            logger.warn("Number of tumor samples isn't equal the number of normal samples");

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(VCFUtils.getHeaderFields(this.getToolkit()));
        headerLines.add(new VCFFormatHeaderLine(SOMATIC_TAG_NAME, 1, VCFHeaderLineType.Float, "Probability that the site is a somatic mutation"));
        headerLines.add(new VCFHeaderLine("source", SOURCE_NAME));
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
    }

    private double log10pNonRefInSamples(final VariantContext vc, final Set<String> samples) {
        return log10pSumInSamples(vc, samples, false);
    }

    private double log10pRefInSamples(final VariantContext vc, final Set<String> samples) {
        return log10pSumInSamples(vc, samples, true);
    }

    private double log10pSumInSamples(final VariantContext vc, final Set<String> samples, boolean calcRefP) {
        double log10p = 0;

        for ( final String sample : samples ) {
            Genotype g = vc.getGenotype(sample);
            if ( g.isNoCall() ) {
                log10p += 0;
            } else {
                double[] gLikelihoods = MathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
                double log10pNonRefSample = Math.log10(calcRefP ? gLikelihoods[0] : 1 - gLikelihoods[0]);
                log10p += log10pNonRefSample;
            }
        }

        return log10p;
    }

    private double calcLog10pSomatic(final VariantContext vc) {
        // walk over tumors, and calculate pNonRef
        double log10pNonRefInTumors = log10pNonRefInSamples(vc, tumorSamples);
        double log10pRefInNormals = log10pRefInSamples(vc, normalSamples);
        double log10SomaticPrior = MathUtils.phredScaleToLog10Probability(somaticPriorQ);
        double log10Somatic = log10SomaticPrior + log10pNonRefInTumors - log10pRefInNormals;
        return log10Somatic;
    }

    /**
     * For each variant in the file, determine the phasing for the child and replace the child's genotype with the trio's genotype
     *
     * @param tracker  the reference meta-data tracker
     * @param ref      the reference context
     * @param context  the alignment context
     * @return null
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            for ( final VariantContext vc : tracker.getValues(variantCollection.variants, context.getLocation()) ) {
                double log10pSomatic = calcLog10pSomatic(vc);

                // write in the somatic status probability
                Map<String, Object> attrs = new HashMap<String, Object>(); // vc.getAttributes());
                attrs.put(SOMATIC_TAG_NAME, MathUtils.log10ProbabilityToPhredScale(log10pSomatic));
                VariantContext newvc = VariantContext.modifyAttributes(vc, attrs);

                vcfWriter.add(newvc);
            }

            return null;
        }

        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Integer reduceInit() {
        return null;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }
}
