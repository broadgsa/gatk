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

package org.broadinstitute.sting.playground.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.Allele;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.*;
import org.broad.tribble.vcf.VCFConstants;
import org.broad.tribble.dbsnp.DbSNPFeature;

import java.io.PrintStream;
import java.util.*;


public class UnifiedGenotyperEngine {

    public static final String TRIGGER_TRACK_NAME = "trigger";
    public static final String LOW_QUAL_FILTER_NAME = "LowQual";

    // the unified argument collection
    private UnifiedArgumentCollection UAC = null;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    // the model used for calculating genotypes
    private ThreadLocal<GenotypeLikelihoodsCalculationModel> glcm = new ThreadLocal<GenotypeLikelihoodsCalculationModel>();

    // the model used for calculating p(non-ref)
    private ThreadLocal<AlleleFrequencyCalculationModel> afcm = new ThreadLocal<AlleleFrequencyCalculationModel>();

    // because the allele frequency priors are constant for a given i, we cache the results to avoid having to recompute everything
    private double[] log10AlleleFrequencyPriors;

    // the allele frequency likelihoods (allocated once as an optimization)
    private ThreadLocal<double[]> log10AlleleFrequencyPosteriors = new ThreadLocal<double[]>();

    // the priors object
    private GenotypePriors genotypePriors;

    // the various loggers and writers
    private Logger logger = null;
    private PrintStream verboseWriter = null;

    // number of chromosomes (2 * samples) in input
    private int N;

    // the standard filter to use for calls below the confidence threshold but above the emit threshold
    private static final Set<String> filter = new HashSet<String>(1);


    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC) {
        // get the number of samples
        // if we're supposed to assume a single sample, do so
        int numSamples;
        if ( UAC.ASSUME_SINGLE_SAMPLE != null )
            numSamples = 1;
        else
            numSamples = SampleUtils.getSAMFileSamples(toolkit.getSAMFileHeader()).size();
        initialize(UAC, null, null, null, numSamples);
    }

    public UnifiedGenotyperEngine(UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, int numSamples) {
        initialize(UAC, logger, verboseWriter, engine, numSamples);
    }

    private void initialize(UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, int numSamples) {
        this.UAC = UAC;
        this.logger = logger;
        this.verboseWriter = verboseWriter;
        this.annotationEngine = engine;

        N = 2 * numSamples;
        log10AlleleFrequencyPriors = new double[N+1];
        computeAlleleFrequencyPriors(N);
        genotypePriors = createGenotypePriors(UAC);

        filter.add(LOW_QUAL_FILTER_NAME);
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public VariantCallContext runGenotyper(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {

        // initialize the GenotypeCalculationModel for this thread if that hasn't been done yet
        if ( glcm.get() == null ) {
            glcm.set(getGenotypeLikelihoodsCalculationObject(logger, UAC));
            log10AlleleFrequencyPosteriors.set(new double[N+1]);
            afcm.set(getAlleleFrequencyCalculationObject(N, logger, verboseWriter, UAC));
        }

        BadReadPileupFilter badReadPileupFilter = new BadReadPileupFilter(refContext);
        Map<String, StratifiedAlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext, badReadPileupFilter);
        if ( stratifiedContexts == null )
            return null;

        Map<String, BiallelicGenotypeLikelihoods> GLs = new HashMap<String, BiallelicGenotypeLikelihoods>();
        Allele refAllele = glcm.get().getLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE, genotypePriors, GLs);

        // estimate our confidence in a reference call and return
        if ( GLs.size() == 0 )
            return estimateReferenceConfidence(stratifiedContexts, genotypePriors.getHeterozygosity(), false);

        // reset the optimization value and determine the p(AF>0)
        // TODO: get rid of this optimization, it is wrong!
        afcm.get().setMinAlleleFrequencyToTest(0);

        // 'zero' out the AFs (so that we don't have to worry if not all samples have reads at this position)
        clearAFarray(log10AlleleFrequencyPosteriors.get());
        afcm.get().getLog10PNonRef(tracker, refContext, GLs, log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.get());

        // find the most likely frequency
        int bestAFguess = MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors.get());

        // calculate p(f>0)
        double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get());
        double sum = 0.0;
        for (int i = 1; i <= N; i++)
            sum += normalizedPosteriors[i];
        double PofF = Math.min(sum, 1.0); // deal with precision errors

        double phredScaledConfidence;
        if ( bestAFguess != 0 ) {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(normalizedPosteriors[0]);
            if ( Double.isInfinite(phredScaledConfidence) )
                // todo - verify this is OK
                phredScaledConfidence = (double)QualityUtils.MAX_QUAL_SCORE;
                //phredScaledConfidence = -10.0 * log10AlleleFrequencyPosteriors.get()[0];
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

        // did we trigger on the provided track?
        boolean atTriggerTrack = tracker.getReferenceMetaData(TRIGGER_TRACK_NAME, false).size() > 0;

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( !UAC.ALL_BASES_MODE && !passesEmitThreshold(phredScaledConfidence, bestAFguess, atTriggerTrack) ) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            return estimateReferenceConfidence(stratifiedContexts, genotypePriors.getHeterozygosity(), true);
        }

        // create the genotypes
        Map<String, Genotype> genotypes = afcm.get().assignGenotypes(stratifiedContexts, GLs, log10AlleleFrequencyPosteriors.get(), bestAFguess);

        // next, get the variant context data (alleles, attributes, etc.)
        HashSet<Allele> alleles = new HashSet<Allele>();
        alleles.add(refAllele);
        for ( Genotype g : genotypes.values() )
            alleles.addAll(g.getAlleles());

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        HashMap<String, Object> attributes = new HashMap<String, Object>();

        DbSNPFeature dbsnp = getDbSNP(tracker);
        if ( dbsnp != null )
            attributes.put(VariantContext.ID_KEY, dbsnp.getRsID());

        // if the site was downsampled, record that fact
        if ( rawContext.hasPileupBeenDownsampled() )
            attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);


        if ( !UAC.NO_SLOD ) {

            // TODO: implement me

            //Map<String, BiallelicGenotypeLikelihoods> forwardGLs = glcm.get().getLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.FORWARD, genotypePriors);
            //Map<String, BiallelicGenotypeLikelihoods> reverseGLs = glcm.get().getLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.REVERSE, genotypePriors);






        }

        GenomeLoc loc = refContext.getLocus();

        // todo - temp fix until we can deal with extended events properly
        //VariantContext vc = new VariantContext("UG_call", loc.getContig(), loc.getStart(), loc.getStop(), alleles, genotypes, phredScaledConfidence/10.0, passesCallThreshold(phredScaledConfidence, atTriggerTrack) ? null : filter, attributes);
        VariantContext vc = new VariantContext("UG_call", loc.getContig(), loc.getStart(),
                Math.max(loc.getStart()+refAllele.length()-1,loc.getStart()),
                alleles, genotypes, phredScaledConfidence/10.0, passesCallThreshold(phredScaledConfidence, atTriggerTrack) ? null : filter, attributes);


        if ( annotationEngine != null ) {
            // first off, we want to use the *unfiltered* context for the annotations
            ReadBackedPileup pileup = null;
            if (rawContext.hasExtendedEventPileup())
                pileup = rawContext.getExtendedEventPileup();
            else if (rawContext.hasBasePileup())
                pileup = rawContext.getBasePileup();

            stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(pileup, UAC.ASSUME_SINGLE_SAMPLE);

            Collection<VariantContext> variantContexts = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vc);
            vc = variantContexts.iterator().next(); //We know the collection will always have exactly 1 element.
        }

        // print out stats if we have a writer
        if ( verboseWriter != null )
            printVerboseData(vc, PofF, phredScaledConfidence, normalizedPosteriors);

        VariantCallContext call = new VariantCallContext(vc, passesCallThreshold(phredScaledConfidence, atTriggerTrack));
        call.setRefBase(refContext.getBase());
        return call;
    }

    private static boolean isValidDeletionFraction(double d) {
        return ( d >= 0.0 && d <= 1.0 );
    }

    private Map<String, StratifiedAlignmentContext> getFilteredAndStratifiedContexts(UnifiedArgumentCollection UAC, ReferenceContext refContext, AlignmentContext rawContext, BadReadPileupFilter badReadPileupFilter) {
        Map<String, StratifiedAlignmentContext> stratifiedContexts = null;

        if ( UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.DINDEL && rawContext.hasExtendedEventPileup() ) {

            ReadBackedExtendedEventPileup rawPileup = rawContext.getExtendedEventPileup();

            // filter the context based on min mapping quality
            ReadBackedExtendedEventPileup pileup = rawPileup.getMappingFilteredPileup(UAC.MIN_MAPPING_QUALTY_SCORE);

            // filter the context based on bad mates and mismatch rate
            pileup = pileup.getFilteredPileup(badReadPileupFilter);

            // don't call when there is no coverage
            if ( pileup.size() == 0 && !UAC.ALL_BASES_MODE )
                return null;

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(pileup, UAC.ASSUME_SINGLE_SAMPLE);

        } else if ( UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.SNP && !rawContext.hasExtendedEventPileup() ) {

            byte ref = refContext.getBase();
            if ( !BaseUtils.isRegularBase(ref) )
                return null;

            ReadBackedPileup rawPileup = rawContext.getBasePileup();

            // filter the context based on min base and mapping qualities
            ReadBackedPileup pileup = rawPileup.getBaseAndMappingFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE, UAC.MIN_MAPPING_QUALTY_SCORE);

            // filter the context based on bad mates and mismatch rate
            pileup = pileup.getFilteredPileup(badReadPileupFilter);

            // don't call when there is no coverage
            if ( pileup.size() == 0 && !UAC.ALL_BASES_MODE )
                return null;

            // are there too many deletions in the pileup?
            if ( isValidDeletionFraction(UAC.MAX_DELETION_FRACTION) &&
                 (double)pileup.getNumberOfDeletions() / (double)pileup.size() > UAC.MAX_DELETION_FRACTION )
                return null;

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(pileup, UAC.ASSUME_SINGLE_SAMPLE);
        }

        return stratifiedContexts;
    }

    protected static void clearAFarray(double[] AFs) {
        for ( int i = 0; i < AFs.length; i++ )
            AFs[i] = AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED;
    }

    private VariantCallContext estimateReferenceConfidence(Map<String, StratifiedAlignmentContext> contexts, double theta, boolean ignoreCoveredSamples) {

        // TODO: implement me

        double P_of_ref = 1.0;

        // use the AF=0 prob if it's calculated
        //if ( ignoreCoveredSamples )
        //    P_of_ref = 1.0 - PofFs[BaseUtils.simpleBaseToBaseIndex(bestAlternateAllele)];

        // for each sample that we haven't examined yet
        //for ( String sample : samples ) {
        //    boolean isCovered = contexts.containsKey(sample);
        //    if ( ignoreCoveredSamples && isCovered )
        //        continue;

        P_of_ref = 0.5;
        //    int depth = isCovered ? contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup().size() : 0;
        //    P_of_ref *= 1.0 - (theta / 2.0) * MathUtils.binomialProbability(0, depth, 0.5);
        //}

        return new VariantCallContext(QualityUtils.phredScaleErrorRate(1.0 - P_of_ref) >= UAC.STANDARD_CONFIDENCE_FOR_CALLING);
    }

    protected void printVerboseData(VariantContext vc, double PofF, double phredScaledConfidence, double[] normalizedPosteriors) {
        for (int i = 0; i <= N; i++) {
            StringBuilder AFline = new StringBuilder("AFINFO\t");
            AFline.append(vc.getChr()).append(":").append(vc.getStart()).append("\t");
            AFline.append(vc.getReference()).append("\t");
            if ( vc.isBiallelic() )
                AFline.append(vc.getAlternateAllele(0)).append("\t");
            else
                AFline.append("N/A\t");
            AFline.append(i + "/" + N + "\t");
            AFline.append(String.format("%.2f\t", ((float)i)/N));
            AFline.append(String.format("%.8f\t", log10AlleleFrequencyPriors[i]));
            if ( log10AlleleFrequencyPosteriors.get()[i] == AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED)
                AFline.append("0.0\t");                
            else
                AFline.append(String.format("%.8f\t", log10AlleleFrequencyPosteriors.get()[i]));
            AFline.append(String.format("%.8f\t", normalizedPosteriors[i]));
            verboseWriter.println(AFline.toString());
        }

        verboseWriter.println("P(f>0) = " + PofF);
        verboseWriter.println("Qscore = " + phredScaledConfidence);
        verboseWriter.println();
    }

    /**
     * Filters low quality reads out of the pileup.
     */
    private class BadReadPileupFilter implements PileupElementFilter {
        private ReferenceContext refContext;

        public BadReadPileupFilter(ReferenceContext ref) {
            // create the +/-20bp window
            GenomeLoc window = GenomeLocParser.createGenomeLoc(ref.getLocus().getContig(), Math.max(ref.getWindow().getStart(), ref.getLocus().getStart()-20), Math.min(ref.getWindow().getStop(), ref.getLocus().getStart()+20));
            byte[] bases = new byte[41];
            System.arraycopy(ref.getBases(), (int)Math.max(0, window.getStart()-ref.getWindow().getStart()), bases, 0, (int)window.size());
            refContext = new ReferenceContext(ref.getLocus(), window, bases);
        }

        public boolean allow(PileupElement pileupElement) {
              return  ((UAC.USE_BADLY_MATED_READS || !BadMateFilter.hasBadMate(pileupElement.getRead())) &&
                     AlignmentUtils.mismatchesInRefWindow(pileupElement, refContext, true) <= UAC.MAX_MISMATCHES );
        }
    }

    /**
     * @param tracker   rod data
     *
     * @return the dbsnp rod if there is one at this position
     */
    protected static DbSNPFeature getDbSNP(RefMetaDataTracker tracker) {
        return DbSNPHelper.getFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
    }

    protected boolean passesEmitThreshold(double conf, int bestAFguess, boolean atTriggerTrack) {
        return (atTriggerTrack ?
                (conf >= Math.min(UAC.TRIGGER_CONFIDENCE_FOR_CALLING, UAC.TRIGGER_CONFIDENCE_FOR_EMITTING)) :
                ((UAC.GENOTYPE_MODE || bestAFguess != 0) && conf >= Math.min(UAC.STANDARD_CONFIDENCE_FOR_CALLING, UAC.STANDARD_CONFIDENCE_FOR_EMITTING)));
    }

    protected boolean passesCallThreshold(double conf, boolean atTriggerTrack) {
        return (atTriggerTrack ?
                (conf >= UAC.TRIGGER_CONFIDENCE_FOR_CALLING) :
                (conf >= UAC.STANDARD_CONFIDENCE_FOR_CALLING));
    }

    protected void computeAlleleFrequencyPriors(int N) {
        // calculate the allele frequency priors for 1-N
        double sum = 0.0;
        for (int i = 1; i <= N; i++) {
            double value = UAC.heterozygosity / (double)i;
            log10AlleleFrequencyPriors[i] = Math.log10(value);
            sum += value;
        }

        // null frequency for AF=0 is (1 - sum(all other frequencies))
        log10AlleleFrequencyPriors[0] = Math.log10(1.0 - sum);
    }

    private static GenotypePriors createGenotypePriors(UnifiedArgumentCollection UAC) {
        GenotypePriors priors;
        switch ( UAC.GLmodel ) {
            case SNP:
                // use flat priors for GLs
                priors = new DiploidSNPGenotypePriors();
                break;
            case DINDEL:
                // create flat priors for Indels, actual priors will depend on event length to be genotyped
                priors = new DiploidIndelGenotypePriors();
                break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + UAC.GLmodel);
        }

        return priors;
    }

    private static GenotypeLikelihoodsCalculationModel getGenotypeLikelihoodsCalculationObject(Logger logger, UnifiedArgumentCollection UAC) {
        GenotypeLikelihoodsCalculationModel glcm;
        switch ( UAC.GLmodel ) {
            case SNP:
                glcm = new SNPGenotypeLikelihoodsCalculationModel(UAC, logger);
                break;
            case DINDEL:
                glcm = new DindelGenotypeLikelihoodsCalculationModel(UAC, logger);
                break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + UAC.GLmodel);
        }

        return glcm;
    }

    private static AlleleFrequencyCalculationModel getAlleleFrequencyCalculationObject(int N, Logger logger, PrintStream verboseWriter, UnifiedArgumentCollection UAC) {
        AlleleFrequencyCalculationModel afcm;
        switch ( UAC.AFmodel ) {
            case EXACT:
                afcm = new ExactAFCalculationModel(N, logger, verboseWriter);
                break;
            case GRID_SEARCH:
                afcm = new GridSearchAFEstimation(N, logger, verboseWriter);
                break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + UAC.GLmodel);
        }

        return afcm;
    }
}