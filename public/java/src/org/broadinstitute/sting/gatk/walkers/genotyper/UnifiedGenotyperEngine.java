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
        /** produces calls only at variant sites */
        EMIT_VARIANTS_ONLY,
        /** produces calls at variant sites and confident reference sites */
        EMIT_ALL_CONFIDENT_SITES,
        /** produces calls at any callable site regardless of confidence; this argument is intended for point
         * mutations (SNPs) only and while some indel calls may be produced they are by no means comprehensive */
        EMIT_ALL_SITES
    }

    protected static final List<Allele> NO_CALL_ALLELES = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
    protected static final double SUM_GL_THRESH_NOCALL = -0.001; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.

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
    private final double[][] log10AlleleFrequencyPriorsSNPs;
    private final double[][] log10AlleleFrequencyPriorsIndels;

    // the allele frequency likelihoods (allocated once as an optimization)
    private ThreadLocal<double[][]> log10AlleleFrequencyLikelihoods = new ThreadLocal<double[][]>();
    private ThreadLocal<double[][]> log10AlleleFrequencyPosteriors = new ThreadLocal<double[][]>();

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

    private final GenomeLocParser genomeLocParser;
    private final boolean BAQEnabledOnCMDLine;



    // ---------------------------------------------------------------------------------------------------------
    //
    // Public interface functions
    //
    // ---------------------------------------------------------------------------------------------------------
    @Requires({"toolkit != null", "UAC != null"})
    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC) {
        this(toolkit, UAC, Logger.getLogger(UnifiedGenotyperEngine.class), null, null, SampleUtils.getSAMFileSamples(toolkit.getSAMFileHeader()));
    }

    @Requires({"toolkit != null", "UAC != null", "logger != null", "samples != null && samples.size() > 0"})
    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, Set<String> samples) {
        this.BAQEnabledOnCMDLine = toolkit.getArguments().BAQMode != BAQ.CalculationMode.OFF;
        genomeLocParser = toolkit.getGenomeLocParser();
        this.samples = new TreeSet<String>(samples);
        // note that, because we cap the base quality by the mapping quality, minMQ cannot be less than minBQ
        this.UAC = UAC.clone();

        this.logger = logger;
        this.verboseWriter = verboseWriter;
        this.annotationEngine = engine;

        N = 2 * this.samples.size();
        log10AlleleFrequencyPriorsSNPs = new double[UAC.MAX_ALTERNATE_ALLELES][N+1];
        log10AlleleFrequencyPriorsIndels = new double[UAC.MAX_ALTERNATE_ALLELES][N+1];
        computeAlleleFrequencyPriors(N, log10AlleleFrequencyPriorsSNPs, UAC.heterozygosity);
        computeAlleleFrequencyPriors(N, log10AlleleFrequencyPriorsIndels, UAC.INDEL_HETEROZYGOSITY);
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

        return glcm.get().get(model).getLikelihoods(tracker, refContext, stratifiedContexts, type, getGenotypePriors(model), alternateAlleleToUse, useBAQedPileup && BAQEnabledOnCMDLine);
    }

    private VariantCallContext generateEmptyContext(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, AlignmentContext rawContext) {
        VariantContext vc;
        if ( UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            VariantContext vcInput = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, rawContext.getLocation(), false, logger, UAC.alleles);
            if ( vcInput == null )
                return null;
            vc = new VariantContextBuilder(vcInput).source("UG_call").noID().referenceBaseForIndel(ref.getBase()).make();
        } else {
            // deal with bad/non-standard reference bases
            if ( !Allele.acceptableAlleleBases(new byte[]{ref.getBase()}) )
                return null;

            Set<Allele> alleles = new HashSet<Allele>();
            alleles.add(Allele.create(ref.getBase(), true));
            vc = new VariantContextBuilder("UG_call", ref.getLocus().getContig(), ref.getLocus().getStart(), ref.getLocus().getStart(), alleles).make();
        }
        
        if ( annotationEngine != null ) {
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            ReadBackedPileup pileup = null;
            if (rawContext.hasExtendedEventPileup())
                pileup = rawContext.getExtendedEventPileup();
            else if (rawContext.hasBasePileup())
                pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);

            vc = annotationEngine.annotateContext(tracker, ref, stratifiedContexts, vc);
        }

        return new VariantCallContext(vc, false);
    }

    public VariantCallContext calculateGenotypes(VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {
        return calculateGenotypes(null, null, null, null, vc, model);
    }

    public VariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {

        boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;

        // initialize the data for this thread if that hasn't been done yet
        if ( afcm.get() == null ) {
            log10AlleleFrequencyLikelihoods.set(new double[UAC.MAX_ALTERNATE_ALLELES][N+1]);
            log10AlleleFrequencyPosteriors.set(new double[UAC.MAX_ALTERNATE_ALLELES][N+1]);
            afcm.set(getAlleleFrequencyCalculationObject(N, logger, verboseWriter, UAC));
        }

        // don't try to genotype too many alternate alleles
        if ( vc.getAlternateAlleles().size() > UAC.MAX_ALTERNATE_ALLELES ) {
            logger.warn("the Unified Genotyper is currently set to genotype at most " + UAC.MAX_ALTERNATE_ALLELES + " alternate alleles in a given context, but the context at " + vc.getChr() + ":" + vc.getStart() + " has " + vc.getAlternateAlleles().size() + " alternate alleles; see the --max_alternate_alleles argument");
            return null;
        }

        // estimate our confidence in a reference call and return
        if ( vc.getNSamples() == 0 ) {
            if ( limitedContext )
                return null;
            return (UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES ?
                    estimateReferenceConfidence(vc, stratifiedContexts, getGenotypePriors(model).getHeterozygosity(), false, 1.0) :
                    generateEmptyContext(tracker, refContext, stratifiedContexts, rawContext));
        }

        // 'zero' out the AFs (so that we don't have to worry if not all samples have reads at this position)
        clearAFarray(log10AlleleFrequencyLikelihoods.get());
        clearAFarray(log10AlleleFrequencyPosteriors.get());
        afcm.get().getLog10PNonRef(vc.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyLikelihoods.get(), log10AlleleFrequencyPosteriors.get());

        // is the most likely frequency conformation AC=0 for all alternate alleles?
        boolean bestGuessIsRef = true;

        // which alternate allele has the highest MLE AC?
        int indexOfHighestAlt = -1;
        int alleleCountOfHighestAlt = -1;

        // determine which alternate alleles have AF>0
        boolean[] altAllelesToUse = new boolean[vc.getAlternateAlleles().size()];
        for ( int i = 0; i < vc.getAlternateAlleles().size(); i++ ) {
            int indexOfBestAC = MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors.get()[i]);

            // if the most likely AC is not 0, then this is a good alternate allele to use
            if ( indexOfBestAC != 0 ) {
                altAllelesToUse[i] = true;
                bestGuessIsRef = false;
            }
            // if in GENOTYPE_GIVEN_ALLELES mode, we still want to allow the use of a poor allele
            else if ( UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
                altAllelesToUse[i] = true;
            }

            // keep track of the "best" alternate allele to use
            if ( indexOfBestAC > alleleCountOfHighestAlt) {
                alleleCountOfHighestAlt = indexOfBestAC;
                indexOfHighestAlt = i;
            }
        }

        // calculate p(f>0)
        // TODO -- right now we just calculate it for the alt allele with highest AF, but the likelihoods need to be combined correctly over all AFs
        double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get()[indexOfHighestAlt]);
        double sum = 0.0;
        for (int i = 1; i <= N; i++)
            sum += normalizedPosteriors[i];
        double PofF = Math.min(sum, 1.0); // deal with precision errors

        double phredScaledConfidence;
        if ( !bestGuessIsRef || UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(normalizedPosteriors[0]);
            if ( Double.isInfinite(phredScaledConfidence) )
                phredScaledConfidence = -10.0 * log10AlleleFrequencyPosteriors.get()[0][0];
        } else {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofF);
            if ( Double.isInfinite(phredScaledConfidence) ) {
                sum = 0.0;
                for (int i = 1; i <= N; i++) {
                    if ( log10AlleleFrequencyPosteriors.get()[0][i] == AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED )
                        break;
                    sum += log10AlleleFrequencyPosteriors.get()[0][i];
                }
                phredScaledConfidence = (MathUtils.compareDoubles(sum, 0.0) == 0 ? 0 : -10.0 * sum);
            }
        }

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES && !passesEmitThreshold(phredScaledConfidence, bestGuessIsRef) ) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            return limitedContext ? null : estimateReferenceConfidence(vc, stratifiedContexts, getGenotypePriors(model).getHeterozygosity(), true, 1.0 - PofF);
        }

        // strip out any alleles that aren't going to be used in the VariantContext
        List<Allele> myAlleles;
        if ( UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY ) {
            myAlleles = new ArrayList<Allele>(vc.getAlleles().size());
            myAlleles.add(vc.getReference());
            for ( int i = 0; i < vc.getAlternateAlleles().size(); i++ ) {
                if ( altAllelesToUse[i] )
                    myAlleles.add(vc.getAlternateAllele(i));
            }
        } else {
            // use all of the alleles if we are given them by the user
            myAlleles = vc.getAlleles();
        }

        // start constructing the resulting VC
        GenomeLoc loc = genomeLocParser.createGenomeLoc(vc);
        VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), loc.getStop(), myAlleles);
        builder.log10PError(phredScaledConfidence/-10.0);
        if ( ! passesCallThreshold(phredScaledConfidence) )
            builder.filters(filter);
        if ( limitedContext ) {
            builder.referenceBaseForIndel(vc.getReferenceBaseForIndel());
        } else {
            builder.referenceBaseForIndel(refContext.getBase());
        }

        // create the genotypes
        GenotypesContext genotypes = assignGenotypes(vc, altAllelesToUse);

        // print out stats if we have a writer
        if ( verboseWriter != null && !limitedContext )
            printVerboseData(refContext.getLocus().toString(), vc, PofF, phredScaledConfidence, normalizedPosteriors, model);

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        HashMap<String, Object> attributes = new HashMap<String, Object>();

        // if the site was downsampled, record that fact
        if ( !limitedContext && rawContext.hasPileupBeenDownsampled() )
            attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);

        if ( UAC.COMPUTE_SLOD && !limitedContext && !bestGuessIsRef ) {
            //final boolean DEBUG_SLOD = false;

            // the overall lod
            VariantContext vcOverall = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, vc.getAlternateAllele(0), false, model);
            clearAFarray(log10AlleleFrequencyLikelihoods.get());
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(vcOverall.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyLikelihoods.get(), log10AlleleFrequencyPosteriors.get());
            //double overallLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double overallLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get()[0], 1);
            //if ( DEBUG_SLOD ) System.out.println("overallLog10PofF=" + overallLog10PofF);

            // the forward lod
            VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.FORWARD, vc.getAlternateAllele(0), false, model);
            clearAFarray(log10AlleleFrequencyLikelihoods.get());
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(vcForward.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyLikelihoods.get(), log10AlleleFrequencyPosteriors.get());
            //double[] normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
            double forwardLog10PofNull = log10AlleleFrequencyPosteriors.get()[0][0];
            double forwardLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get()[0], 1);
            //if ( DEBUG_SLOD ) System.out.println("forwardLog10PofNull=" + forwardLog10PofNull + ", forwardLog10PofF=" + forwardLog10PofF);

            // the reverse lod
            VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.REVERSE, vc.getAlternateAllele(0), false, model);
            clearAFarray(log10AlleleFrequencyLikelihoods.get());
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(vcReverse.getGenotypes(), vc.getAlleles(), getAlleleFrequencyPriors(model), log10AlleleFrequencyLikelihoods.get(), log10AlleleFrequencyPosteriors.get());
            //normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
            double reverseLog10PofNull = log10AlleleFrequencyPosteriors.get()[0][0];
            double reverseLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get()[0], 1);
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

        // finish constructing the resulting VC
        builder.genotypes(genotypes);
        builder.attributes(attributes);
        VariantContext vcCall = builder.make();

        if ( annotationEngine != null && !limitedContext ) {
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            ReadBackedPileup pileup = null;
            if (rawContext.hasExtendedEventPileup())
                pileup = rawContext.getExtendedEventPileup();
            else if (rawContext.hasBasePileup())
                pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);

            vcCall = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vcCall);
        }

        return new VariantCallContext(vcCall, confidentlyCalled(phredScaledConfidence, PofF));
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
                stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);

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
                stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);
            }
        } else if ( model == GenotypeLikelihoodsCalculationModel.Model.SNP ) {

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(rawContext.getBasePileup());

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

    protected static void clearAFarray(double[][] AFs) {
        for ( int i = 0; i < AFs.length; i++ ) {
            for ( int j = 0; j < AFs[i].length; j++ ) {
                AFs[i][j] = AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED;
            }
        }
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
            if ( log10AlleleFrequencyPosteriors.get()[0][i] == AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED)
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

    protected boolean passesEmitThreshold(double conf, boolean bestGuessIsRef) {
        return (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) && conf >= Math.min(UAC.STANDARD_CONFIDENCE_FOR_CALLING, UAC.STANDARD_CONFIDENCE_FOR_EMITTING);
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

    protected static void computeAlleleFrequencyPriors(final int N, final double[][] priors, final double theta) {

        // the dimension here is the number of alternate alleles; with e.g. 2 alternate alleles the prior will be theta^2 / i
        for (int alleles = 1; alleles <= priors.length; alleles++) {
            double sum = 0.0;

            // for each i
            for (int i = 1; i <= N; i++) {
                double value = Math.pow(theta, alleles) / (double)i;
                priors[alleles-1][i] = Math.log10(value);
                sum += value;
            }

            // null frequency for AF=0 is (1 - sum(all other frequencies))
            priors[alleles-1][0] = Math.log10(1.0 - sum);
        }
    }

    protected double[][] getAlleleFrequencyPriors( final GenotypeLikelihoodsCalculationModel.Model model ) {
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

    /**
     * @param vc            variant context with genotype likelihoods
     * @param allelesToUse  bit vector describing which alternate alleles from the vc are okay to use
     *
     * @return genotypes
     */
    public GenotypesContext assignGenotypes(VariantContext vc,
                                            boolean[] allelesToUse) {

        final GenotypesContext GLs = vc.getGenotypes();

        final List<String> sampleIndices = GLs.getSampleNamesOrderedByName();

        final GenotypesContext calls = GenotypesContext.create();

        for ( int k = GLs.size() - 1; k >= 0; k-- ) {
            final String sample = sampleIndices.get(k);
            final Genotype g = GLs.get(sample);
            if ( !g.hasLikelihoods() )
                continue;

            final double[] likelihoods = g.getLikelihoods().getAsVector();

            // if there is no mass on the likelihoods, then just no-call the sample
            if ( MathUtils.sum(likelihoods) > SUM_GL_THRESH_NOCALL ) {
                calls.add(new Genotype(g.getSampleName(), NO_CALL_ALLELES, Genotype.NO_LOG10_PERROR, null, null, false));
                continue;
            }

            // genotype likelihoods are a linear vector that can be thought of as a row-wise upper triangular matrix of log10Likelihoods.
            // so e.g. with 2 alt alleles the likelihoods are AA,AB,AC,BB,BC,CC and with 3 alt alleles they are AA,AB,AC,AD,BB,BC,BD,CC,CD,DD.

            final int numAltAlleles = allelesToUse.length;

            // start with the assumption that the ideal genotype is homozygous reference
            Allele maxAllele1 = vc.getReference(), maxAllele2 = vc.getReference();
            double maxLikelihoodSeen = likelihoods[0];
            int indexOfMax = 0;

            // keep track of some state
            Allele firstAllele = vc.getReference();
            int subtractor = numAltAlleles + 1;
            int subtractionsMade = 0;

            for ( int i = 1, PLindex = 1; i < likelihoods.length; i++, PLindex++ ) {
                if ( PLindex == subtractor ) {
                    firstAllele = vc.getAlternateAllele(subtractionsMade);
                    PLindex -= subtractor;
                    subtractor--;
                    subtractionsMade++;

                    // we can skip this allele if it's not usable
                    if ( !allelesToUse[subtractionsMade-1] ) {
                        i += subtractor - 1;
                        PLindex += subtractor - 1;
                        continue;
                    }
                }

                // we don't care about the entry if we've already seen better
                if ( likelihoods[i] <= maxLikelihoodSeen )
                    continue;

                // if it's usable then update the alleles
                int alleleIndex = subtractionsMade + PLindex - 1;
                if ( allelesToUse[alleleIndex] ) {
                    maxAllele1 = firstAllele;
                    maxAllele2 = vc.getAlternateAllele(alleleIndex);
                    maxLikelihoodSeen = likelihoods[i];
                    indexOfMax = i;
                }
            }

            ArrayList<Allele> myAlleles = new ArrayList<Allele>();
            myAlleles.add(maxAllele1);
            myAlleles.add(maxAllele2);

            final double qual = GenotypeLikelihoods.getQualFromLikelihoods(indexOfMax, likelihoods);
            calls.add(new Genotype(sample, myAlleles, qual, null, g.getAttributes(), false));
        }

        return calls;
    }
}
