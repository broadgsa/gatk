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

    @Argument(shortName="somaticMinLOD", fullName="somaticMinLOD", required=false, doc="Phred-scaled min probability that a site should be called somatic mutation")
    public byte somaticMinLOD = 1;

    @Output
    protected VCFWriter vcfWriter = null;

    private final String SOMATIC_LOD_TAG_NAME = "SOMATIC_LOD";
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
        headerLines.add(new VCFInfoHeaderLine(VCFConstants.SOMATIC_KEY, 0, VCFHeaderLineType.Flag, "Is this a confidently called somatic mutation"));
        headerLines.add(new VCFFormatHeaderLine(SOMATIC_LOD_TAG_NAME, 1, VCFHeaderLineType.Float, "log10 probability that the site is a somatic mutation"));
        headerLines.add(new VCFHeaderLine("source", SOURCE_NAME));
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
    }

    private double log10pNonRefInSamples(final VariantContext vc, final Set<String> samples) {
        double[] log10ps = log10PLFromSamples(vc, samples, false);
        return MathUtils.log10sumLog10(log10ps); // product of probs => prod in real space
     }

    private double log10pRefInSamples(final VariantContext vc, final Set<String> samples) {
        double[] log10ps = log10PLFromSamples(vc, samples, true);
        return MathUtils.sum(log10ps); // product is sum
    }

    private double[] log10PLFromSamples(final VariantContext vc, final Set<String> samples, boolean calcRefP) {
        double[] log10p = new double[samples.size()];

        int i = 0;
        for ( final String sample : samples ) {
            Genotype g = vc.getGenotype(sample);
            double log10pSample = -1000;
            if ( ! g.isNoCall() ) {
                double[] gLikelihoods = MathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
                log10pSample = Math.log10(calcRefP ? gLikelihoods[0] : 1 - gLikelihoods[0]);
                log10pSample = Double.isInfinite(log10pSample) ? -10000 : log10pSample;
            }
            log10p[i++] = log10pSample;
        }

        return log10p;
    }

    /**
     * P(somatic | D)
     *   = P(somatic) * P(D | somatic)
     *   = P(somatic) * P(D | normals are ref) * P(D | tumors are non-ref)
     *
     * P(! somatic | D)
     *   = P(! somatic) * P(D | ! somatic)
     *   = P(! somatic) *
     *      * (  P(D | normals are non-ref) * P(D | tumors are non-ref) [germline]
     *         + P(D | normals are ref) * P(D | tumors are ref)) [no-variant at all]
     *
     * @param vc
     * @return
     */
    private double calcLog10pSomatic(final VariantContext vc) {
        // walk over tumors
        double log10pNonRefInTumors = log10pNonRefInSamples(vc, tumorSamples);
        double log10pRefInTumors = log10pRefInSamples(vc, tumorSamples);

        // walk over normals
        double log10pNonRefInNormals = log10pNonRefInSamples(vc, normalSamples);
        double log10pRefInNormals = log10pRefInSamples(vc, normalSamples);

        // priors
        double log10pSomaticPrior = MathUtils.phredScaleToLog10Probability(somaticPriorQ);
        double log10pNotSomaticPrior = Math.log10(1 - MathUtils.phredScaleToProbability(somaticPriorQ));

        double log10pNotSomaticGermline = log10pNonRefInNormals + log10pNonRefInTumors;
        double log10pNotSomaticNoVariant = log10pRefInNormals + log10pRefInTumors;

        double log10pNotSomatic = log10pNotSomaticPrior + MathUtils.log10sumLog10(new double[]{log10pNotSomaticGermline, log10pNotSomaticNoVariant});
        double log10pSomatic = log10pSomaticPrior + log10pNonRefInTumors + log10pRefInNormals;
        double lod = log10pSomatic - log10pNotSomatic;

        return Double.isInfinite(lod) ? -10000 : lod;
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
                attrs.put(SOMATIC_LOD_TAG_NAME, log10pSomatic);
                if ( log10pSomatic > somaticMinLOD )
                    attrs.put(VCFConstants.SOMATIC_KEY, true);
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
