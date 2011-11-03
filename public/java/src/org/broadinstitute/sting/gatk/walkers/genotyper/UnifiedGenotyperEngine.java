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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class UnifiedGenotyperEngine {
    public static final String LOW_QUAL_FILTER_NAME = "LowQual";

    public enum OUTPUT_MODE {
        /** the default */
        EMIT_VARIANTS_ONLY,
        /** include confident reference sites */
        EMIT_ALL_CONFIDENT_SITES,
        /** any callable site regardless of confidence */
        EMIT_ALL_SITES
    }

    // the unified argument collection
    private final UnifiedArgumentCollection UAC;
    public UnifiedArgumentCollection getUAC() { return UAC; }

    // the annotation engine
    private final VariantAnnotatorEngine annotationEngine;

    // the model used for calculating genotypes
    private ThreadLocal<Map<GenotypeLikelihoodsCalculationModel.Model, GenotypeLikelihoodsCalculationModel>> glcm = new ThreadLocal<Map<GenotypeLikelihoodsCalculationModel.Model, GenotypeLikelihoodsCalculationModel>>();

    // the model used for calculating p(non-ref)
    private ThreadLocal<AlleleFrequencyCalculationModel> afcm = new ThreadLocal<AlleleFrequencyCalculationModel>();

    // because the allele frequency priors are constant for a given i, we cache the results to avoid having to recompute everything
    private final double[] log10AlleleFrequencyPriorsSNPs;
    private final double[] log10AlleleFrequencyPriorsIndels;

    // the allele frequency likelihoods (allocated once as an optimization)
    private ThreadLocal<double[]> log10AlleleFrequencyPosteriors = new ThreadLocal<double[]>();

    // the priors object
    private final GenotypePriors genotypePriorsSNPs;
    private final GenotypePriors genotypePriorsIndels;

    // samples in input
    private final Set<String> samples;

    // the various loggers and writers
    private final Logger logger;
    private final PrintStream verboseWriter;

    // number of chromosomes (2 * samples) in input
    private final int N;

    // the standard filter to use for calls below the confidence threshold but above the emit threshold
    private static final Set<String> filter = new HashSet<String>(1);

    private final boolean BAQEnabledOnCMDLine;



    // ---------------------------------------------------------------------------------------------------------
    //
    // Public interface functions
    //
    // ---------------------------------------------------------------------------------------------------------
    @Requires({"toolkit != null", "UAC != null"})
    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC) {
        this(toolkit, UAC, Logger.getLogger(UnifiedGenotyperEngine.class), null, null,
                // get the number of samples
                // if we're supposed to assume a single sample, do so
                UAC.ASSUME_SINGLE_SAMPLE != null ?
                        new TreeSet<String>(Arrays.asList(UAC.ASSUME_SINGLE_SAMPLE)) :
                        SampleUtils.getSAMFileSamples(toolkit.getSAMFileHeader()));
    }

    @Requires({"toolkit != null", "UAC != null", "logger != null", "samples != null && samples.size() > 0"})
    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, Set<String> samples) {
        this.BAQEnabledOnCMDLine = toolkit.getArguments().BAQMode != BAQ.CalculationMode.OFF;
        this.samples = new TreeSet<String>(samples);
        // note that, because we cap the base quality by the mapping quality, minMQ cannot be less than minBQ
        this.UAC = UAC.clone();

        this.logger = logger;
        this.verboseWriter = verboseWriter;
        this.annotationEngine = engine;

        N = 2 * this.samples.size();
        log10AlleleFrequencyPriorsSNPs = new double[N+1];
        log10AlleleFrequencyPriorsIndels = new double[N+1];
        computeAlleleFrequencyPriors(N, log10AlleleFrequencyPriorsSNPs, GenotypeLikelihoodsCalculationModel.Model.SNP);
        computeAlleleFrequencyPriors(N, log10AlleleFrequencyPriorsIndels, GenotypeLikelihoodsCalculationModel.Model.INDEL);
        genotypePriorsSNPs = createGenotypePriors(GenotypeLikelihoodsCalculationModel.Model.SNP);
        genotypePriorsIndels = createGenotypePriors(GenotypeLikelihoodsCalculationModel.Model.INDEL);
        
        filter.add(LOW_QUAL_FILTER_NAME);
    }

    /**
     * Compute full calls at a given locus. Entry point for engine calls from the UnifiedGenotyper.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateLikelihoodsAndGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        final GenotypeLikelihoodsCalculationModel.Model model = getCurrentGLModel(tracker, refContext, rawContext );
        if( model == null ) {
            return (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(tracker, refContext, null, rawContext) : null);
        }

        Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext, model);
        if ( stratifiedContexts == null ) {
            return (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(tracker, refContext, stratifiedContexts, rawContext) : null);
        }
        
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, null, true, model);

        if ( vc == null )
            return null;

        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model);
    }

    /**
     * Compute GLs at a given locus. Entry point for engine calls from UGCalcLikelihoods.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantContext object
     */
    public VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        final GenotypeLikelihoodsCalculationModel.Model model = getCurrentGLModel( tracker, refContext, rawContext );
        if( model == null )
            return null;

        Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext, model);
        if ( stratifiedContexts == null )
            return null;

        return calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, null, true, model);
    }

    /**
     * Compute genotypes at a given locus. Entry point for engine calls from UGCallVariants.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @param vc         the GL-annotated variant context
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, VariantContext vc) {
        final GenotypeLikelihoodsCalculationModel.Model model = getCurrentGLModel(tracker, refContext, rawContext );
        if( model == null ) {
            return null;
        }
        Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext, model);
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model);
    }




    // ---------------------------------------------------------------------------------------------------------
    //
    // Private implementation helpers
    //
    // ---------------------------------------------------------------------------------------------------------

    // private method called by both UnifiedGenotyper and UGCalcLikelihoods entry points into the engine
    private VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, Map<String, AlignmentContext> stratifiedContexts, AlignmentContextUtils.ReadOrientation type, Allele alternateAlleleToUse, boolean useBAQedPileup, final GenotypeLikelihoodsCalculationModel.Model model) {

        // initialize the data for this thread if that hasn't been done yet
        if ( glcm.get() == null ) {
            glcm.set(getGenotypeLikelihoodsCalculationObject(logger, UAC));
        }

        Map<String, MultiallelicGenotypeLikelihoods> GLs = new HashMap<String, MultiallelicGenotypeLikelihoods>();

        Allele refAllele = glcm.get().get(model).getLikelihoods(tracker, refContext, stratifiedContexts, type, getGenotypePriors(model), GLs, alternateAlleleToUse, useBAQedPileup && BAQEnabledOnCMDLine);

        if ( refAllele != null )
            return createVariantContextFromLikelihoods(refContext, refAllele, GLs);
        else
            return null;
    }

    private VariantCallContext generateEmptyContext(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, AlignmentContext rawContext) {
        VariantContext vc;
        if ( UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            VariantContext vcInput = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, rawContext.getLocation(), false, logger, UAC.alleles);
            if ( vcInput == null )
                return null;
            vc = new VariantContext("UG_call", vcInput.getChr(), vcInput.getStart(), vcInput.getEnd(), vcInput.getAlleles(), InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null, ref.getBase());

        } else {
            // deal with bad/non-standard reference bases
            if ( !Allele.acceptableAlleleBases(new byte[]{ref.getBase()}) )
                return null;

            Set<Allele> alleles = new HashSet<Allele>();
            alleles.add(Allele.create(ref.getBase(), true));
            vc = new VariantContext("UG_call", ref.getLocus().getContig(), ref.getLocus().getStart(), ref.getLocus().getStart(), alleles);
        }
        
        if ( annotationEngine != null ) {
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            ReadBackedPileup pileup = null;
            if (rawContext.hasExtendedEventPileup())
                pileup = rawContext.getExtendedEventPileup();
            else if (rawContext.hasBasePileup())
                pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup, UAC.ASSUME_SINGLE_SAMPLE);

            vc = annotationEngine.annotateContext(tracker, ref, stratifiedContexts, vc);
        }

        return new VariantCallContext(vc, false);
    }

    private VariantContext createVariantContextFromLikelihoods(ReferenceContext refContext, Allele refAllele, Map<String, MultiallelicGenotypeLikelihoods> GLs) {
        // no-call everyone for now
        List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        Set<Allele> alleles = new LinkedHashSet<Allele>();
        alleles.add(refAllele);
        boolean addedAltAlleles = false;

        HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
        for ( MultiallelicGenotypeLikelihoods GL : GLs.values() ) {
            if ( !addedAltAlleles ) {
                addedAltAlleles = true;
                // ordering important to maintain consistency
                for (Allele a: GL.getAlleles()) {
                    alleles.add(a);
                }
            }

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            //GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
            GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(GL.getLikelihoods());
            attributes.put(VCFConstants.DEPTH_KEY, GL.getDepth());
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, likelihoods);

            genotypes.put(GL.getSample(), new Genotype(GL.getSample(), noCall, Genotype.NO_NEG_LOG_10PERROR, null, attributes, false));
        }

        GenomeLoc loc = refContext.getLocus();
        int endLoc = calculateEndPos(alleles, refAllele, loc);

        return new VariantContext("UG_call",
                loc.getContig(),
                loc.getStart(),
                endLoc,
                alleles,
                genotypes,
                VariantContext.NO_NEG_LOG_10PERROR,
                null,
                null,
                refContext.getBase());
    }

    // private method called by both UnifiedGenotyper and UGCallVariants entry points into the engine
    private VariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {

        // initialize the data for this thread if that hasn't been done yet
        if ( afcm.get() == null ) {
            log10AlleleFrequencyPosteriors.set(new double[N+1]);
            afcm.set(getAlleleFrequencyCalculationObject(N, logger, verboseWriter, UAC));
        }

        // estimate our confidence in a reference call and return
        if ( vc.getNSamples() == 0 )
            return (UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES ?
                    estimateReferenceConfidence(vc, stratifiedContexts, getGenotypePriors(model).getHeterozygosity(), false, 1.0) :
                    generateEmptyContext(tracker, refContext, stratifiedContexts, rawContext));

        // 'zero' out the AFs (so that we don't have to worry if not all samples have reads at this position)
        clearAFarray(log10AlleleFrequencyPosteriors.get());
        afcm.get().getLog10PNonRef(vc.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyPosteriors.get());

        // find the most likely frequency
        int bestAFguess = MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors.get());

        // calculate p(f>0)
        double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get());
        double sum = 0.0;
        for (int i = 1; i <= N; i++)
            sum += normalizedPosteriors[i];
        double PofF = Math.min(sum, 1.0); // deal with precision errors

        double phredScaledConfidence;
        if ( bestAFguess != 0 || UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(normalizedPosteriors[0]);
            if ( Double.isInfinite(phredScaledConfidence) )
                phredScaledConfidence = -10.0 * log10AlleleFrequencyPosteriors.get()[0];
        } else {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofF);
            if ( Double.isInfinite(phredScaledConfidence) ) {
                sum = 0.0;
                for (int i = 1; i <= N; i++) {
                    if ( log10AlleleFrequencyPosteriors.get()[i] == AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED )
                        break;
                    sum += log10AlleleFrequencyPosteriors.get()[i];
                }
                phredScaledConfidence = (MathUtils.compareDoubles(sum, 0.0) == 0 ? 0 : -10.0 * sum);
            }
        }

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES && !passesEmitThreshold(phredScaledConfidence, bestAFguess) ) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            return estimateReferenceConfidence(vc, stratifiedContexts, getGenotypePriors(model).getHeterozygosity(), true, 1.0 - PofF);
        }

        // create the genotypes
        Map<String, Genotype> genotypes = afcm.get().assignGenotypes(vc, log10AlleleFrequencyPosteriors.get(), bestAFguess);

        // print out stats if we have a writer
        if ( verboseWriter != null )
            printVerboseData(refContext.getLocus().toString(), vc, PofF, phredScaledConfidence, normalizedPosteriors, model);

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        HashMap<String, Object> attributes = new HashMap<String, Object>();

        // if the site was downsampled, record that fact
        if ( rawContext.hasPileupBeenDownsampled() )
            attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);


        if ( UAC.COMPUTE_SLOD && bestAFguess != 0 ) {
            //final boolean DEBUG_SLOD = false;

            // the overall lod
            VariantContext vcOverall = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, vc.getAlternateAllele(0), false, model);
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(vcOverall.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyPosteriors.get());
            //double overallLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double overallLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            //if ( DEBUG_SLOD ) System.out.println("overallLog10PofF=" + overallLog10PofF);

            // the forward lod
            VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.FORWARD, vc.getAlternateAllele(0), false, model);
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(vcForward.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyPosteriors.get());
            //double[] normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
            double forwardLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double forwardLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            //if ( DEBUG_SLOD ) System.out.println("forwardLog10PofNull=" + forwardLog10PofNull + ", forwardLog10PofF=" + forwardLog10PofF);

            // the reverse lod
            VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.REVERSE, vc.getAlternateAllele(0), false, model);
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(vcReverse.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyPosteriors.get());
            //normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
            double reverseLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double reverseLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            //if ( DEBUG_SLOD ) System.out.println("reverseLog10PofNull=" + reverseLog10PofNull + ", reverseLog10PofF=" + reverseLog10PofF);

            double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
            double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;
            //if ( DEBUG_SLOD ) System.out.println("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);

            // strand score is max bias between forward and reverse strands
            double strandScore = Math.max(forwardLod, reverseLod);
            // rescale by a factor of 10
            strandScore *= 10.0;
            //logger.debug(String.format("SLOD=%f", strandScore));

            attributes.put("SB", strandScore);
        }

        GenomeLoc loc = refContext.getLocus();

        int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);

        Set<Allele> myAlleles = new HashSet<Allele>(vc.getAlleles());
        // strip out the alternate allele if it's a ref call
        if ( bestAFguess == 0 && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY ) {
            myAlleles = new HashSet<Allele>(1);
            myAlleles.add(vc.getReference());
        }
        VariantContext vcCall = new VariantContext("UG_call", loc.getContig(), loc.getStart(), endLoc,
                myAlleles, genotypes, phredScaledConfidence/10.0, passesCallThreshold(phredScaledConfidence) ? null : filter, attributes, refContext.getBase());

        if ( annotationEngine != null ) {
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            ReadBackedPileup pileup = null;
            if (rawContext.hasExtendedEventPileup())
                pileup = rawContext.getExtendedEventPileup();
            else if (rawContext.hasBasePileup())
                pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup, UAC.ASSUME_SINGLE_SAMPLE);

            vcCall = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vcCall);
        }

        return new VariantCallContext(vcCall, confidentlyCalled(phredScaledConfidence, PofF));
    }

    // A barebones entry point to the exact model when there is no tracker or stratified contexts available -- only GLs
    public VariantCallContext calculateGenotypes(final VariantContext vc, final GenomeLoc loc, final GenotypeLikelihoodsCalculationModel.Model model) {

        // initialize the data for this thread if that hasn't been done yet
        if ( afcm.get() == null ) {
            log10AlleleFrequencyPosteriors.set(new double[N+1]);
            afcm.set(getAlleleFrequencyCalculationObject(N, logger, verboseWriter, UAC));
        }

        // estimate our confidence in a reference call and return
        if ( vc.getNSamples() == 0 )
            return null;

        // 'zero' out the AFs (so that we don't have to worry if not all samples have reads at this position)
        clearAFarray(log10AlleleFrequencyPosteriors.get());
        afcm.get().getLog10PNonRef(vc.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyPosteriors.get());

        // find the most likely frequency
        int bestAFguess = MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors.get());

        // calculate p(f>0)
        double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get());
        double sum = 0.0;
        for (int i = 1; i <= N; i++)
            sum += normalizedPosteriors[i];
        double PofF = Math.min(sum, 1.0); // deal with precision errors

        double phredScaledConfidence;
        if ( bestAFguess != 0 || UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(normalizedPosteriors[0]);
            if ( Double.isInfinite(phredScaledConfidence) )
                phredScaledConfidence = -10.0 * log10AlleleFrequencyPosteriors.get()[0];
        } else {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofF);
            if ( Double.isInfinite(phredScaledConfidence) ) {
                sum = 0.0;
                for (int i = 1; i <= N; i++) {
                    if ( log10AlleleFrequencyPosteriors.get()[i] == AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED )
                        break;
                    sum += log10AlleleFrequencyPosteriors.get()[i];
                }
                phredScaledConfidence = (MathUtils.compareDoubles(sum, 0.0) == 0 ? 0 : -10.0 * sum);
            }
        }

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES && !passesEmitThreshold(phredScaledConfidence, bestAFguess) ) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            return null;
        }

        // create the genotypes
        Map<String, Genotype> genotypes = afcm.get().assignGenotypes(vc, log10AlleleFrequencyPosteriors.get(), bestAFguess);

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        HashMap<String, Object> attributes = new HashMap<String, Object>();

        int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);

        Set<Allele> myAlleles = new HashSet<Allele>(vc.getAlleles());
        // strip out the alternate allele if it's a ref call
        if ( bestAFguess == 0 && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY ) {
            myAlleles = new HashSet<Allele>(1);
            myAlleles.add(vc.getReference());
        }
        VariantContext vcCall = new VariantContext("UG_call", loc.getContig(), loc.getStart(), endLoc,
                myAlleles, genotypes, phredScaledConfidence/10.0, passesCallThreshold(phredScaledConfidence) ? null : filter, attributes, vc.getReferenceBaseForIndel());

        return new VariantCallContext(vcCall, confidentlyCalled(phredScaledConfidence, PofF));
    }

    private int calculateEndPos(Collection<Allele> alleles, Allele refAllele, GenomeLoc loc) {
        // TODO - temp fix until we can deal with extended events properly
        // for indels, stop location is one more than ref allele length
        boolean isSNP = true, hasNullAltAllele = false;
        for (Allele a : alleles){
            if (a.length() != 1) {
                isSNP = false;
                break;
            }
        }
        for (Allele a : alleles){
            if (a.isNull()) {
                hasNullAltAllele = true;
                break;
            }
        }
        // standard deletion: ref allele length = del length. endLoc = startLoc + refAllele.length(), alt allele = null
        // standard insertion: ref allele length = 0, endLos = startLoc
        // mixed: want end loc = start Loc for case {A*,AT,T} but say  {ATG*,A,T} : want then end loc = start loc + refAllele.length
        // So, in general, end loc = startLoc + refAllele.length, except in complex substitutions where it's one less
        //
        // todo - this is unnecessarily complicated and is so just because of Tribble's arbitrary vc conventions, should be cleaner/simpler,
        // the whole vc processing infrastructure seems too brittle and riddled with special case handling


        int endLoc = loc.getStart();
        if ( !isSNP) {
            endLoc += refAllele.length();
            if(!hasNullAltAllele)
                endLoc--;

        }

        return endLoc;
    }

    private Map<String, AlignmentContext> getFilteredAndStratifiedContexts(UnifiedArgumentCollection UAC, ReferenceContext refContext, AlignmentContext rawContext, final GenotypeLikelihoodsCalculationModel.Model model) {

        Map<String, AlignmentContext> stratifiedContexts = null;

        if ( !BaseUtils.isRegularBase( refContext.getBase() ) )
            return null;

        if ( model == GenotypeLikelihoodsCalculationModel.Model.INDEL ) {

            if (UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES) {
                // regular pileup in this case
                ReadBackedPileup pileup = rawContext.getBasePileup() .getMappingFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE);

                // don't call when there is no coverage
                if ( pileup.getNumberOfElements() == 0 && UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES  )
                    return null;

                // stratify the AlignmentContext and cut by sample
                stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup, UAC.ASSUME_SINGLE_SAMPLE);

            } else {

                // todo - tmp will get rid of extended events so this wont be needed
                if (!rawContext.hasExtendedEventPileup())
                    return null;
                ReadBackedExtendedEventPileup rawPileup = rawContext.getExtendedEventPileup();

                // filter the context based on min mapping quality
                ReadBackedExtendedEventPileup pileup = rawPileup.getMappingFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE);

                // don't call when there is no coverage
                if ( pileup.getNumberOfElements() == 0 && UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES  )
                    return null;

                // stratify the AlignmentContext and cut by sample
                stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup, UAC.ASSUME_SINGLE_SAMPLE);
            }
        } else if ( model == GenotypeLikelihoodsCalculationModel.Model.SNP ) {

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(rawContext.getBasePileup(), UAC.ASSUME_SINGLE_SAMPLE);

            if( !(UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES) ) {
                int numDeletions = 0;
                for( final PileupElement p : rawContext.getBasePileup() ) {
                    if( p.isDeletion() ) { numDeletions++; }
                }
                if( ((double) numDeletions) / ((double) rawContext.getBasePileup().getNumberOfElements()) > UAC.MAX_DELETION_FRACTION ) {
                    return null;
                }
            }
        }

        return stratifiedContexts;
    }

    protected static void clearAFarray(double[] AFs) {
        for ( int i = 0; i < AFs.length; i++ )
            AFs[i] = AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED;
    }

    private final static double[] binomialProbabilityDepthCache = new double[10000];
    static {
        for ( int i = 1; i < binomialProbabilityDepthCache.length; i++ ) {
            binomialProbabilityDepthCache[i] = MathUtils.binomialProbability(0, i, 0.5);
        }
    }

    private final double getRefBinomialProb(final int depth) {
        if ( depth < binomialProbabilityDepthCache.length )
            return binomialProbabilityDepthCache[depth];
        else
            return MathUtils.binomialProbability(0, depth, 0.5);
    }


    private VariantCallContext estimateReferenceConfidence(VariantContext vc, Map<String, AlignmentContext> contexts, double theta, boolean ignoreCoveredSamples, double initialPofRef) {
        if ( contexts == null )
            return null;

        double P_of_ref = initialPofRef;

        // for each sample that we haven't examined yet
        for ( String sample : samples ) {
            boolean isCovered = contexts.containsKey(sample);
            if ( ignoreCoveredSamples && isCovered )
                continue;


            int depth = 0;

            if (isCovered) {
                AlignmentContext context =  contexts.get(sample);
                if (context.hasBasePileup())
                    depth = context.getBasePileup().depthOfCoverage();
                else if (context.hasExtendedEventPileup())
                    depth = context.getExtendedEventPileup().depthOfCoverage();
            }

            P_of_ref *= 1.0 - (theta / 2.0) * getRefBinomialProb(depth);
        }

        return new VariantCallContext(vc, QualityUtils.phredScaleErrorRate(1.0 - P_of_ref) >= UAC.STANDARD_CONFIDENCE_FOR_CALLING, false);
    }

    protected void printVerboseData(String pos, VariantContext vc, double PofF, double phredScaledConfidence, double[] normalizedPosteriors, final GenotypeLikelihoodsCalculationModel.Model model) {
        Allele refAllele = null, altAllele = null;
        for ( Allele allele : vc.getAlleles() ) {
            if ( allele.isReference() )
                refAllele = allele;
            else
                altAllele = allele;
        }

        for (int i = 0; i <= N; i++) {
            StringBuilder AFline = new StringBuilder("AFINFO\t");
            AFline.append(pos);
            AFline.append("\t");
            AFline.append(refAllele);
            AFline.append("\t");
            if ( altAllele != null )
                AFline.append(altAllele);
            else
                AFline.append("N/A");
            AFline.append("\t");
            AFline.append(i + "/" + N + "\t");
            AFline.append(String.format("%.2f\t", ((float)i)/N));
            AFline.append(String.format("%.8f\t", getAlleleFrequencyPriors(model)[i]));
            if ( log10AlleleFrequencyPosteriors.get()[i] == AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED)
                AFline.append("0.00000000\t");                
            else
                AFline.append(String.format("%.8f\t", log10AlleleFrequencyPosteriors.get()[i]));
            AFline.append(String.format("%.8f\t", normalizedPosteriors[i]));
            verboseWriter.println(AFline.toString());
        }

        verboseWriter.println("P(f>0) = " + PofF);
        verboseWriter.println("Qscore = " + phredScaledConfidence);
        verboseWriter.println();
    }

    protected boolean passesEmitThreshold(double conf, int bestAFguess) {
        return (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES || bestAFguess != 0) && conf >= Math.min(UAC.STANDARD_CONFIDENCE_FOR_CALLING, UAC.STANDARD_CONFIDENCE_FOR_EMITTING);
    }

    protected boolean passesCallThreshold(double conf) {
        return conf >= UAC.STANDARD_CONFIDENCE_FOR_CALLING;
    }

    protected boolean confidentlyCalled(double conf, double PofF) {
        return conf >= UAC.STANDARD_CONFIDENCE_FOR_CALLING ||
                (UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES && QualityUtils.phredScaleErrorRate(PofF) >= UAC.STANDARD_CONFIDENCE_FOR_CALLING);
    }

    // decide whether we are currently processing SNPs, indels, or neither
    private GenotypeLikelihoodsCalculationModel.Model getCurrentGLModel(final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                                        final AlignmentContext rawContext ) {
        if (rawContext.hasExtendedEventPileup() ) {
            // todo - remove this code
            if ((UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BOTH || UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.INDEL) &&
                   (UAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) )
                return GenotypeLikelihoodsCalculationModel.Model.INDEL;
        }
        else {
            // no extended event pileup
            // if we're genotyping given alleles and we have a requested SNP at this position, do SNP
            if (UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES) {
                VariantContext vcInput = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, refContext, rawContext.getLocation(), false, logger, UAC.alleles);
                if (vcInput == null)
                    return null;

                // todo - no support to genotype MNP's yet
                if  (vcInput.isMNP())
                    return null;

                if (vcInput.isSNP())  {
                    if (( UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BOTH || UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.SNP))
                        return GenotypeLikelihoodsCalculationModel.Model.SNP;
                    else
                        // ignore SNP's if user chose INDEL mode
                        return null;
                }
                else if ((vcInput.isIndel() || vcInput.isMixed()) && (UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BOTH || UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.INDEL))
                    return GenotypeLikelihoodsCalculationModel.Model.INDEL;
            }
            else {
                // todo - this assumes SNP's take priority when BOTH is selected, should do a smarter way once extended events are removed
                if( UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BOTH || UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.SNP)
                    return GenotypeLikelihoodsCalculationModel.Model.SNP;
                else if (UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.INDEL)
                    return GenotypeLikelihoodsCalculationModel.Model.INDEL;
            }
        }
        return null;
    }

    protected void computeAlleleFrequencyPriors(int N, final double[] priors, final GenotypeLikelihoodsCalculationModel.Model model) {
        // calculate the allele frequency priors for 1-N
        double sum = 0.0;
        double heterozygosity;

        if (model == GenotypeLikelihoodsCalculationModel.Model.INDEL)
            heterozygosity = UAC.INDEL_HETEROZYGOSITY;
        else
            heterozygosity = UAC.heterozygosity;
        
        for (int i = 1; i <= N; i++) {
            double value = heterozygosity / (double)i;
            priors[i] = Math.log10(value);
            sum += value;
        }

        // null frequency for AF=0 is (1 - sum(all other frequencies))
        priors[0] = Math.log10(1.0 - sum);
    }

    protected double[] getAlleleFrequencyPriors( final GenotypeLikelihoodsCalculationModel.Model model ) {
        switch( model ) {
            case SNP:
                return log10AlleleFrequencyPriorsSNPs;
            case INDEL:
                return log10AlleleFrequencyPriorsIndels;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);
        }
    }

    private static GenotypePriors createGenotypePriors( final GenotypeLikelihoodsCalculationModel.Model model ) {
        GenotypePriors priors;
        switch ( model ) {
            case SNP:
                // use flat priors for GLs
                priors = new DiploidSNPGenotypePriors();
                break;
            case INDEL:
                // create flat priors for Indels, actual priors will depend on event length to be genotyped
                priors = new DiploidIndelGenotypePriors();
                break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);
        }
        return priors;
    }

    protected GenotypePriors getGenotypePriors( final GenotypeLikelihoodsCalculationModel.Model model ) {
        switch( model ) {
            case SNP:
                return genotypePriorsSNPs;
            case INDEL:
                return genotypePriorsIndels;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);
        }
    }

    private static Map<GenotypeLikelihoodsCalculationModel.Model, GenotypeLikelihoodsCalculationModel> getGenotypeLikelihoodsCalculationObject(Logger logger, UnifiedArgumentCollection UAC) {
        Map<GenotypeLikelihoodsCalculationModel.Model, GenotypeLikelihoodsCalculationModel> glcm = new HashMap<GenotypeLikelihoodsCalculationModel.Model, GenotypeLikelihoodsCalculationModel>();
        glcm.put(GenotypeLikelihoodsCalculationModel.Model.SNP, new SNPGenotypeLikelihoodsCalculationModel(UAC, logger));
        glcm.put(GenotypeLikelihoodsCalculationModel.Model.INDEL, new IndelGenotypeLikelihoodsCalculationModel(UAC, logger));
        return glcm;
    }

    private static AlleleFrequencyCalculationModel getAlleleFrequencyCalculationObject(int N, Logger logger, PrintStream verboseWriter, UnifiedArgumentCollection UAC) {
        AlleleFrequencyCalculationModel afcm;
        switch ( UAC.AFmodel ) {
            case EXACT:
                afcm = new ExactAFCalculationModel(UAC, N, logger, verboseWriter);
                break;
            case GRID_SEARCH:
                afcm = new GridSearchAFEstimation(UAC, N, logger, verboseWriter);
                break;
            default: throw new IllegalArgumentException("Unexpected AlleleFrequencyCalculationModel " + UAC.AFmodel);
        }

        return afcm;
    }

    public static VariantContext getVCFromAllelesRod(RefMetaDataTracker tracker, ReferenceContext ref, GenomeLoc loc, boolean requireSNP, Logger logger, final RodBinding<VariantContext> allelesBinding) {
        if ( tracker == null || ref == null || logger == null )
            throw new ReviewedStingException("Bad arguments: tracker=" + tracker + " ref=" + ref + " logger=" + logger);
        VariantContext vc = null;

        // search for usable record
        for( final VariantContext vc_input : tracker.getValues(allelesBinding, loc) ) {
            if ( vc_input != null && ! vc_input.isFiltered() && (! requireSNP || vc_input.isSNP() )) {
                if ( vc == null ) {
                    vc = vc_input;
                } else {
                    logger.warn("Multiple valid VCF records detected in the alleles input file at site " + ref.getLocus() + ", only considering the first record");
                }
            }
        }

        return vc;
    }
}
