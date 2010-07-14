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

package org.broadinstitute.sting.playground.gatk.walkers;

import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.MutableVariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.gatk.walkers.varianteval.MendelianViolationEvaluator;

import java.util.*;
import java.io.File;

/**
 * Implements an (STILL BEING TESTED) algorithm for calling SNPs in trios
 */

@By(DataSource.REFERENCE)
//@Requires(value={DataSource.REFERENCE, DataSource.REFERENCE_BASES, DataSource.READS},referenceMetaData={@RMD(name="sites",type= VariationRod.class)})
@Allows({DataSource.READS, DataSource.REFERENCE})
//, @RMD(name="parent1",type= VariationRod.class), @RMD(name="parent2",type= VariationRod.class)})
public class TrioGenotyperWalker extends RefWalker<VariantContext, Integer>{
    @Argument(shortName="mom", doc="", required=true)
    protected String mom;

    @Argument(shortName="dad", doc="", required=true)
    protected String dad;

    @Argument(shortName="kid", doc="", required=true)
    protected String kid;

    @Argument(shortName="log10PriorOfDeNovoOfTrueVariant", doc="", required=false)
    double LOG10_MENDEL_VIOLATION_PRIOR = -5;   // 30 in 3B bases

    @Argument(shortName = "varout", doc = "File to which variants should be written", required = true)
    public String vcfOutputFile = null;

    @ArgumentCollection
    private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    UnifiedGenotyperEngine UGEngine = null;
    private List<String> FAMILY_MEMBERS;
    private VCFWriter writer = null;

    public void initialize() {
        UGEngine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, null);
        // initialize the header
        FAMILY_MEMBERS = Arrays.asList(mom, dad, kid);

        // initialize the writer
        writer = new VCFWriter(new File(vcfOutputFile));
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariantContext vc = tracker.getVariantContext(ref, "variants", EnumSet.of(VariantContext.Type.SNP), context.getLocation(), true);
        
        if ( vc != null && vc.isPolymorphic() ) {
            if ( ! vc.hasGenotypes(FAMILY_MEMBERS) )
                throw new StingException("variants file does not contain genotypes for everyone in family: " + FAMILY_MEMBERS);

            VariantCallContext call = UGEngine.runGenotyper(tracker, ref, context);

            // is call ever be null?
            vc = annotateTrioPrior(vc, call.vc);

            return vc;
        } else {
            return null;
        }
    }

    private VariantContext annotateTrioPrior(VariantContext vcIn, VariantContext call) {
        Genotype momG = call.getGenotype(mom);
        Genotype dadG = call.getGenotype(dad);
        Genotype kidG = call.getGenotype(kid);

        double log10POfGenotype = Double.MIN_VALUE;
        if ( MendelianViolationEvaluator.isViolation(call, momG, dadG, kidG) ) {
            Allele R = call.getReference();
            Allele A = call.getAlternateAllele(0);

            List<List<Allele>> possibleGenotypes = Arrays.asList(Arrays.asList(R,R), Arrays.asList(R,A), Arrays.asList(A,A));

            double[] L = new double[3 * 3 * 3];
            int i = 0, bestIndex = 0;
            double log10LOfBestGenotypes = Integer.MIN_VALUE;
            for ( List<Allele> momPG : possibleGenotypes ) {
                for ( List<Allele> dadPG : possibleGenotypes ) {
                    for ( List<Allele> kidPG : possibleGenotypes ) {
                        double log10LOfG = genotypeL(momPG, momG) + genotypeL(dadPG, dadG) + genotypeL(kidPG, kidG);
                        boolean isViolation = MendelianViolationEvaluator.isViolation(call, momPG, dadPG, kidPG);
                        double log10prior = isViolation ? LOG10_MENDEL_VIOLATION_PRIOR : 0;
                        L[i] = log10LOfG + log10prior;

                        if ( log10LOfG > log10LOfBestGenotypes ) {
                            bestIndex = i;
                            log10LOfBestGenotypes = log10LOfG;
                        }
                        logger.debug(String.format("%10s %10s => %10s : %b\t%.2f\t\t%.2f\t\t%3.2f", momPG, dadPG, kidPG, isViolation, log10LOfG, log10prior, L[i]));
                        i++;
                    }
                }
            }

            double[] posteriors = MathUtils.normalizeFromLog10(L, true);
            log10POfGenotype = posteriors[bestIndex];
        }
        //log10POfViolation = Math.min(log10POfViolation, 0);

        double Q = QualityUtils.phredScaleCorrectRate(Math.pow(10, log10POfGenotype));
        logger.debug(String.format("log10 P of best genotype log10 post = %.2f, Q = %.2f", log10POfGenotype, Q));
        MutableVariantContext mvc = new MutableVariantContext(vcIn);
        mvc.putAttribute("MVQ", Q);
        return new VariantContext(mvc);
    }

    /**
     * Isolate the rest of the walker from the code to get genotype likelihood values for allele A/B in genotypeCall
     * @param alleles
     * @param genotypeCall
     * @return
     */
    private double genotypeL( List<Allele> alleles, Genotype genotypeCall ) {
        String postTriplet = (String)genotypeCall.getAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY);
        if ( postTriplet == null )
            throw new StingException("BUG: TrioGenotyperWalker expected genotype likelihood triplets " + VCFConstants.GENOTYPE_LIKELIHOODS_KEY);

        // calculate the offset -- AA => 0, AB => 1, BB => 2
        int i = 0;
        for ( Allele a : alleles )
            i += a.isNonReference() ? 1 : 0;

        // convert the corresponding GL field to a double
        String log10LStrings[] = postTriplet.split(",");
        return Double.valueOf(log10LStrings[i]);
    }

    public Integer reduceInit() { return 0; }
    public Integer reduce(VariantContext vc, Integer a) {
        if ( vc != null ) {
            if ( a == 0 )
                writer.writeHeader(VariantContextAdaptors.createVCFHeader(null, vc));

            writer.add(vc, new byte[]{(byte)'.'});
            a++;
        }

        return a;
    }

    public void onTraversalDone(Integer result) {} // Don't print the reduce result
}