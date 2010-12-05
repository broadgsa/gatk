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

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.GenotypeLikelihoods;
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
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.sam.GATKSAMRecordFilter;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broad.tribble.vcf.VCFConstants;

import java.io.PrintStream;
import java.util.*;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.reference.IndexedFastaSequenceFile;


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

    // samples in input
    private Set<String> samples = new TreeSet<String>();

    // the various loggers and writers
    private Logger logger = null;
    private PrintStream verboseWriter = null;

    // fasta reference reader to supplement the edges of the reference sequence for long reads
    private IndexedFastaSequenceFile referenceReader;

    // number of chromosomes (2 * samples) in input
    private int N;

    // the standard filter to use for calls below the confidence threshold but above the emit threshold
    private static final Set<String> filter = new HashSet<String>(1);

    private static final int MISMATCH_WINDOW_SIZE = 20;

    
    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC) {
        // get the number of samples
        // if we're supposed to assume a single sample, do so
        samples = new TreeSet<String>();
        if ( UAC.ASSUME_SINGLE_SAMPLE != null )
            samples.add(UAC.ASSUME_SINGLE_SAMPLE);
        else
            samples = SampleUtils.getSAMFileSamples(toolkit.getSAMFileHeader());

        initialize(toolkit, UAC, null, null, null, samples.size());
    }

    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, Set<String> samples) {
        this.samples = new TreeSet<String>(samples);
        initialize(toolkit, UAC, logger, verboseWriter, engine, samples.size());
    }

    private void initialize(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, int numSamples) {
        // note that, because we cap the base quality by the mapping quality, minMQ cannot be less than minBQ
        this.UAC = UAC.clone();
        this.UAC.MIN_MAPPING_QUALTY_SCORE = Math.max(UAC.MIN_MAPPING_QUALTY_SCORE, UAC.MIN_BASE_QUALTY_SCORE);

        this.logger = logger;
        this.verboseWriter = verboseWriter;
        this.annotationEngine = engine;

        N = 2 * numSamples;
        log10AlleleFrequencyPriors = new double[N+1];
        computeAlleleFrequencyPriors(N);
        genotypePriors = createGenotypePriors(UAC);

        filter.add(LOW_QUAL_FILTER_NAME);

        referenceReader = new CachingIndexedFastaSequenceFile(toolkit.getArguments().referenceFile);
    }

    /**
     * Compute full calls at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateLikelihoodsAndGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        Map<String, StratifiedAlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext);
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE, null);
        if ( vc == null )
            return null;

        VariantCallContext vcc = calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc);
        vcc.vc = GLsToPLs(vcc.vc);
        return vcc;
    }

    /**
     * Compute GLs at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @param alternateAlleleToUse the alternate allele to use, null if not set
     * @return the VariantContext object
     */
    public VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Allele alternateAlleleToUse) {
        Map<String, StratifiedAlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext);
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE, alternateAlleleToUse);
        return GLsToPLs(vc);
    }

    private VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType type, Allele alternateAlleleToUse) {
        if ( stratifiedContexts == null )
            return null;

        // initialize the data for this thread if that hasn't been done yet
        if ( glcm.get() == null ) {
            glcm.set(getGenotypeLikelihoodsCalculationObject(logger, UAC));
        }

        Map<String, BiallelicGenotypeLikelihoods> GLs = new HashMap<String, BiallelicGenotypeLikelihoods>();
        Allele refAllele = glcm.get().getLikelihoods(tracker, refContext, stratifiedContexts, type, genotypePriors, GLs, alternateAlleleToUse);

        if (refAllele != null)
            return createVariantContextFromLikelihoods(refContext, refAllele, GLs);
        else
            return null;
    }

    private VariantContext createVariantContextFromLikelihoods(ReferenceContext refContext, Allele refAllele, Map<String, BiallelicGenotypeLikelihoods> GLs) {
        // no-call everyone for now
        List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        Set<Allele> alleles = new HashSet<Allele>();
        alleles.add(refAllele);
        boolean addedAltAllele = false;

        HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
        for ( BiallelicGenotypeLikelihoods GL : GLs.values() ) {
            if ( !addedAltAllele ) {
                addedAltAllele = true;
                alleles.add(GL.getAlleleA());
                alleles.add(GL.getAlleleB());
            }

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods(), VCFConstants.GENOTYPE_LIKELIHOODS_KEY);
            attributes.put(VCFConstants.DEPTH_KEY, GL.getDepth());
            attributes.put(likelihoods.getKey(), likelihoods.getAsString());

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
                null);
    }

    private static VariantContext GLsToPLs(VariantContext vc) {
        if ( vc == null )
            return null;

        HashMap<String, Genotype> calls = new HashMap<String, Genotype>();
        for ( Map.Entry<String, Genotype> genotype : vc.getGenotypes().entrySet() ) {
            HashMap<String, Object> attributes = new HashMap<String, Object>(genotype.getValue().getAttributes());
            attributes.remove(VCFConstants.GENOTYPE_LIKELIHOODS_KEY);
            attributes.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.GLsToPLs(genotype.getValue().getLikelihoods().getAsVector()));
            calls.put(genotype.getKey(), Genotype.modifyAttributes(genotype.getValue(), attributes));
        }
        return VariantContext.modifyGenotypes(vc, calls);
    }

    /**
     * Compute genotypes at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @param vc         the GL-annotated variant context
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, VariantContext vc) {
        Map<String, StratifiedAlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext);
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc);
    }

    private VariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        // initialize the data for this thread if that hasn't been done yet
        if ( afcm.get() == null ) {
            log10AlleleFrequencyPosteriors.set(new double[N+1]);
            afcm.set(getAlleleFrequencyCalculationObject(N, logger, verboseWriter, UAC));
        }

        // estimate our confidence in a reference call and return
        if ( vc.getNSamples() == 0 )
            return estimateReferenceConfidence(stratifiedContexts, genotypePriors.getHeterozygosity(), false, 1.0);

        // 'zero' out the AFs (so that we don't have to worry if not all samples have reads at this position)
        clearAFarray(log10AlleleFrequencyPosteriors.get());
        afcm.get().getLog10PNonRef(tracker, refContext, vc.getGenotypes(), log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.get());

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

        // did we trigger on the provided track?
        boolean atTriggerTrack = tracker.getReferenceMetaData(TRIGGER_TRACK_NAME, false).size() > 0;

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( !UAC.ALL_BASES_MODE && !passesEmitThreshold(phredScaledConfidence, bestAFguess, atTriggerTrack) ) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            return estimateReferenceConfidence(stratifiedContexts, genotypePriors.getHeterozygosity(), true, 1.0 - PofF);
        }

        // create the genotypes
        Map<String, Genotype> genotypes = afcm.get().assignGenotypes(vc, log10AlleleFrequencyPosteriors.get(), bestAFguess);

        // print out stats if we have a writer
        if ( verboseWriter != null )
            printVerboseData(refContext.getLocus().toString(), vc, PofF, phredScaledConfidence, normalizedPosteriors);

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        HashMap<String, Object> attributes = new HashMap<String, Object>();

        String rsID = DbSNPHelper.rsIDOfFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
        if ( rsID != null )
            attributes.put(VariantContext.ID_KEY, rsID);

        // if the site was downsampled, record that fact
        if ( rawContext.hasPileupBeenDownsampled() )
            attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);


        if ( !UAC.NO_SLOD && bestAFguess != 0 ) {
            final boolean DEBUG_SLOD = false;

            // the overall lod
            //double overallLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double overallLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            if ( DEBUG_SLOD ) System.out.println("overallLog10PofF=" + overallLog10PofF);

            // the forward lod
            VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.FORWARD, vc.getAlternateAllele(0));
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(tracker, refContext, vcForward.getGenotypes(), log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.get());
            //double[] normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
            double forwardLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double forwardLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            if ( DEBUG_SLOD ) System.out.println("forwardLog10PofNull=" + forwardLog10PofNull + ", forwardLog10PofF=" + forwardLog10PofF);

            // the reverse lod
            VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.REVERSE, vc.getAlternateAllele(0));
            clearAFarray(log10AlleleFrequencyPosteriors.get());
            afcm.get().getLog10PNonRef(tracker, refContext, vcReverse.getGenotypes(), log10AlleleFrequencyPriors, log10AlleleFrequencyPosteriors.get());
            //normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
            double reverseLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
            double reverseLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors.get(), 1);
            if ( DEBUG_SLOD ) System.out.println("reverseLog10PofNull=" + reverseLog10PofNull + ", reverseLog10PofF=" + reverseLog10PofF);

            double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
            double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;
            if ( DEBUG_SLOD ) System.out.println("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);

            // strand score is max bias between forward and reverse strands
            double strandScore = Math.max(forwardLod, reverseLod);
            // rescale by a factor of 10
            strandScore *= 10.0;
            //logger.debug(String.format("SLOD=%f", strandScore));

            attributes.put("SB", Double.valueOf(strandScore));
        }

        GenomeLoc loc = refContext.getLocus();

        int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);

        VariantContext vcCall = new VariantContext("UG_call", loc.getContig(), loc.getStart(), endLoc,
                vc.getAlleles(), genotypes, phredScaledConfidence/10.0, passesCallThreshold(phredScaledConfidence, atTriggerTrack) ? null : filter, attributes);

        if ( annotationEngine != null ) {
            // first off, we want to use the *unfiltered* context for the annotations
            ReadBackedPileup pileup = null;
            if (rawContext.hasExtendedEventPileup())
                pileup = rawContext.getExtendedEventPileup();
            else if (rawContext.hasBasePileup())
                pileup = rawContext.getBasePileup();
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySampleName(pileup, UAC.ASSUME_SINGLE_SAMPLE);

            Collection<VariantContext> variantContexts = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vcCall);
            vcCall = variantContexts.iterator().next(); // we know the collection will always have exactly 1 element.
        }

        VariantCallContext call = new VariantCallContext(vcCall, passesCallThreshold(phredScaledConfidence, atTriggerTrack));
        call.setRefBase(refContext.getBase());
        return call;
    }

    private int calculateEndPos(Set<Allele> alleles, Allele refAllele, GenomeLoc loc) {
        // TODO - temp fix until we can deal with extended events properly
        // for indels, stop location is one more than ref allele length
        boolean isSNP = true;
        for (Allele a : alleles){
            if (a.getBaseString().length() != 1) {
                isSNP = false;
                break;
            }
        }

        int endLoc = loc.getStart();
        if ( !isSNP )
            endLoc += refAllele.length();

        return endLoc;
    }

    private static boolean isValidDeletionFraction(double d) {
        return ( d >= 0.0 && d <= 1.0 );
    }

    private Map<String, StratifiedAlignmentContext> getFilteredAndStratifiedContexts(UnifiedArgumentCollection UAC, ReferenceContext refContext, AlignmentContext rawContext) {
        BadBaseFilter badReadPileupFilter = new BadBaseFilter(refContext, UAC);

        Map<String, StratifiedAlignmentContext> stratifiedContexts = null;

        if ( UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.DINDEL && rawContext.hasExtendedEventPileup() ) {

            ReadBackedExtendedEventPileup rawPileup = rawContext.getExtendedEventPileup();

            // filter the context based on min mapping quality
            ReadBackedExtendedEventPileup pileup = rawPileup.getMappingFilteredPileup(UAC.MIN_MAPPING_QUALTY_SCORE);

            // don't call when there is no coverage
            if ( pileup.size() == 0 && !UAC.ALL_BASES_MODE )
                return null;

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySampleName(pileup, UAC.ASSUME_SINGLE_SAMPLE);

        } else if ( UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.SNP && !rawContext.hasExtendedEventPileup() ) {

            byte ref = refContext.getBase();
            if ( !BaseUtils.isRegularBase(ref) )
                return null;

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySampleName(rawContext.getBasePileup(), UAC.ASSUME_SINGLE_SAMPLE);

            // filter the reads (and test for bad pileups)
            if ( !filterPileup(stratifiedContexts, badReadPileupFilter) )
                return null;
        }

        return stratifiedContexts;
    }

    protected static void clearAFarray(double[] AFs) {
        for ( int i = 0; i < AFs.length; i++ )
            AFs[i] = AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED;
    }

    private VariantCallContext estimateReferenceConfidence(Map<String, StratifiedAlignmentContext> contexts, double theta, boolean ignoreCoveredSamples, double initialPofRef) {
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
                AlignmentContext context =  contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

                if (context.hasBasePileup())
                    depth = context.getBasePileup().size();
                else if (context.hasExtendedEventPileup())
                    depth = context.getExtendedEventPileup().size();
            }

            P_of_ref *= 1.0 - (theta / 2.0) * MathUtils.binomialProbability(0, depth, 0.5);
        }

        return new VariantCallContext(QualityUtils.phredScaleErrorRate(1.0 - P_of_ref) >= UAC.STANDARD_CONFIDENCE_FOR_CALLING);
    }

    protected void printVerboseData(String pos, VariantContext vc, double PofF, double phredScaledConfidence, double[] normalizedPosteriors) {
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
            AFline.append(String.format("%.8f\t", log10AlleleFrequencyPriors[i]));
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

    private boolean filterPileup(Map<String, StratifiedAlignmentContext> stratifiedContexts, BadBaseFilter badBaseFilter) {
        int numDeletions = 0, pileupSize = 0;

        for ( StratifiedAlignmentContext context : stratifiedContexts.values() ) {
            ReadBackedPileup pileup = context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            for ( PileupElement p : pileup ) {
                final SAMRecord read = p.getRead();

                if ( p.isDeletion() ) {
                    // if it's a good read, count it
                    if ( read.getMappingQuality() >= UAC.MIN_MAPPING_QUALTY_SCORE &&
                         (UAC.USE_BADLY_MATED_READS || !BadMateFilter.hasBadMate(read)) )
                        numDeletions++;
                } else {
                    if ( !(read instanceof GATKSAMRecord) )
                        throw new ReviewedStingException("The UnifiedGenotyper currently expects GATKSAMRecords, but instead saw a " + read.getClass());
                    GATKSAMRecord GATKrecord = (GATKSAMRecord)read;
                    GATKrecord.setGoodBases(badBaseFilter, true);
                    if ( GATKrecord.isGoodBase(p.getOffset()) )
                        pileupSize++;
                }
            }
        }

        // now, test for bad pileups

        // in all_bases mode, it doesn't matter
        if ( UAC.ALL_BASES_MODE )
            return true;

        // is there no coverage?
        if ( pileupSize == 0 )
            return false;

        // are there too many deletions in the pileup?
        if ( isValidDeletionFraction(UAC.MAX_DELETION_FRACTION) &&
                (double)numDeletions / (double)(pileupSize + numDeletions) > UAC.MAX_DELETION_FRACTION )
            return false;

        return true;
    }

    /**
     * Filters low quality bases out of the SAMRecords.
     */
    private class BadBaseFilter implements GATKSAMRecordFilter {
        private ReferenceContext refContext;
        private final UnifiedArgumentCollection UAC;

        public BadBaseFilter(ReferenceContext refContext, UnifiedArgumentCollection UAC) {
            this.refContext = refContext;
            this.UAC = UAC;
        }

        public BitSet getGoodBases(final GATKSAMRecord record) {
            // all bits are set to false by default
            BitSet bitset = new BitSet(record.getReadLength());

            // if the mapping quality is too low or the mate is bad, we can just zero out the whole read and continue
            if ( record.getMappingQuality() < UAC.MIN_MAPPING_QUALTY_SCORE ||
                 (!UAC.USE_BADLY_MATED_READS && BadMateFilter.hasBadMate(record)) ) {
                return bitset;
            }

            byte[] quals = record.getBaseQualities();
            for (int i = 0; i < quals.length; i++) {
                if ( quals[i] >= UAC.MIN_BASE_QUALTY_SCORE )
                    bitset.set(i);
            }

            // if a read is too long for the reference context, extend the context (being sure not to extend past the end of the chromosome)
            if ( record.getAlignmentEnd() > refContext.getWindow().getStop() ) {
                GenomeLoc window = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getWindow().getStart(), Math.min(record.getAlignmentEnd(), referenceReader.getSequenceDictionary().getSequence(refContext.getLocus().getContig()).getSequenceLength()));
                byte[] bases = referenceReader.getSubsequenceAt(window.getContig(), window.getStart(), window.getStop()).getBases();
                StringUtil.toUpperCase(bases);
                refContext = new ReferenceContext(refContext.getGenomeLocParser(),refContext.getLocus(), window, bases);
            }            

            BitSet mismatches = AlignmentUtils.mismatchesInRefWindow(record, refContext, UAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE);
            if ( mismatches != null )
                bitset.and(mismatches);

            return bitset;
        }
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