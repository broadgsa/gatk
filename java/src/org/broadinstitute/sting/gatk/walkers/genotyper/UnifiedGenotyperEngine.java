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
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.pileup.*;
import org.broad.tribble.vcf.VCFConstants;

import java.io.PrintStream;
import java.util.*;


public class UnifiedGenotyperEngine {

    public static final String TRIGGER_TRACK_NAME = "trigger";
    public static final String LOW_QUAL_FILTER_NAME = "LowQual";

    protected HashMap<String, String> dbAnnotations = new HashMap<String, String>();

    // the unified argument collection
    protected UnifiedArgumentCollection UAC = null;

    // the annotation engine
    protected VariantAnnotatorEngine annotationEngine;

    // the model used for calculating genotypes
    protected ThreadLocal<GenotypeCalculationModel> gcm = new ThreadLocal<GenotypeCalculationModel>();

    // the various loggers and writers
    protected Logger logger = null;
    protected VCFWriter vcfWriter = null;
    protected PrintStream verboseWriter = null;

    // samples in input
    protected Set<String> samples = new HashSet<String>();


    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC) {
        initialize(toolkit, UAC, null, null, null, null);
    }

    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, VCFWriter genotypeWriter, PrintStream verboseWriter, VariantAnnotatorEngine engine) {
        initialize(toolkit, UAC, logger, genotypeWriter, verboseWriter, engine);

    }

    private void initialize(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, VCFWriter genotypeWriter, PrintStream verboseWriter, VariantAnnotatorEngine engine) {
        this.UAC = UAC;
        this.logger = logger;
        this.vcfWriter = genotypeWriter;
        this.verboseWriter = verboseWriter;
        this.annotationEngine = engine;

        // deal with input errors
        if ( toolkit.getArguments().numberOfThreads > 1 && UAC.ASSUME_SINGLE_SAMPLE != null ) {
            // the ASSUME_SINGLE_SAMPLE argument can't be handled (at least for now) while we are multi-threaded because the IO system doesn't know how to get the sample name
            throw new IllegalArgumentException("For technical reasons, the ASSUME_SINGLE_SAMPLE argument cannot be used with multiple threads; please run again without the -nt argument");
        }

        // get all of the unique sample names
        // if we're supposed to assume a single sample, do so
        if ( UAC.ASSUME_SINGLE_SAMPLE != null )
            this.samples.add(UAC.ASSUME_SINGLE_SAMPLE);
        else
            this.samples = SampleUtils.getSAMFileSamples(toolkit.getSAMFileHeader());

        // check to see whether comp rods were included
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            if ( source.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                dbAnnotations.put(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME, VCFConstants.DBSNP_KEY);
            }
            else if ( source.getName().startsWith(VariantAnnotatorEngine.dbPrefix) ) {
                dbAnnotations.put(source.getName(), source.getName().substring(VariantAnnotatorEngine.dbPrefix.length()));
            }
        }
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
        if ( gcm.get() == null ) {
            gcm.set(GenotypeCalculationModelFactory.makeGenotypeCalculation(samples, logger, UAC, verboseWriter));
        }

        byte ref = refContext.getBase();
        if ( !BaseUtils.isRegularBase(ref) )
            return null;

        VariantCallContext call;
        BadReadPileupFilter badReadPileupFilter = new BadReadPileupFilter(refContext);

        if ( rawContext.hasExtendedEventPileup() ) {

            ReadBackedExtendedEventPileup rawPileup = rawContext.getExtendedEventPileup();

            // filter the context based on min mapping quality
            ReadBackedExtendedEventPileup pileup = rawPileup.getMappingFilteredPileup(UAC.MIN_MAPPING_QUALTY_SCORE);

            // filter the context based on bad mates and mismatch rate
            pileup = pileup.getFilteredPileup(badReadPileupFilter);

            // don't call when there is no coverage
            if ( pileup.size() == 0 && !UAC.ALL_BASES_MODE )
                return null;

            // stratify the AlignmentContext and cut by sample
            Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(pileup, UAC.ASSUME_SINGLE_SAMPLE);
            if ( stratifiedContexts == null )
                return null;

            call = gcm.get().callExtendedLocus(tracker, refContext.getBasesAtLocus(-1), rawContext.getLocation(), stratifiedContexts);

        } else {

            ReadBackedPileup rawPileup = rawContext.getBasePileup();

            // filter the context based on min base and mapping qualities
            ReadBackedPileup pileup = rawPileup.getBaseAndMappingFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE, UAC.MIN_MAPPING_QUALTY_SCORE);

            // filter the context based on bad mates and mismatch rate
            pileup = pileup.getFilteredPileup(badReadPileupFilter);

            // don't call when there is no coverage
            if ( pileup.size() == 0 && !UAC.ALL_BASES_MODE )
                return null;

            // are there too many deletions in the pileup?
            if ( UAC.genotypeModel != GenotypeCalculationModel.Model.DINDEL &&
                 isValidDeletionFraction(UAC.MAX_DELETION_FRACTION) &&
                 (double)pileup.getNumberOfDeletions() / (double)pileup.size() > UAC.MAX_DELETION_FRACTION )
                return null;

            // stratify the AlignmentContext and cut by sample
            Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(pileup, UAC.ASSUME_SINGLE_SAMPLE);
            if ( stratifiedContexts == null )
                return null;

            DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, UAC.heterozygosity, DiploidGenotypePriors.PROB_OF_REFERENCE_ERROR);
            call = gcm.get().callLocus(tracker, ref, rawContext.getLocation(), stratifiedContexts, priors);

            // annotate the call, if possible
            if ( call != null && call.vc != null && annotationEngine != null ) {
                // first off, we want to use the *unfiltered* context for the annotations
                stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(rawContext.getBasePileup());

                Collection<VariantContext> variantContexts = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, call.vc);
                call.vc = variantContexts.iterator().next(); //We know the collection will always have exactly 1 element.
            }
        }

        if ( call != null && call.vc != null ) {
            call.setRefBase(ref);

            // if the site was downsampled, record that fact
            if ( rawContext.hasPileupBeenDownsampled() ) {
                Map<String, Object> attrs = new HashMap<String, Object>(call.vc.getAttributes());
                attrs.put(VCFConstants.DOWNSAMPLED_KEY, true);
                VariantContextUtils.modifyAttributes(call.vc, attrs);
            }
        }
        return call;
    }

    private static boolean isValidDeletionFraction(double d) {
        return ( d >= 0.0 && d <= 1.0 );
    }

    /**
     * Filters low quality reads out of the pileup.
     */
    private class BadReadPileupFilter implements PileupElementFilter {
        private ReferenceContext refContext;

        public BadReadPileupFilter(ReferenceContext refContext) { this.refContext = refContext; }

        public boolean allow(PileupElement pileupElement) {
            return  ((UAC.USE_BADLY_MATED_READS || !BadMateFilter.hasBadMate(pileupElement.getRead())) &&
                     AlignmentUtils.mismatchesInRefWindow(pileupElement, refContext, true) <= UAC.MAX_MISMATCHES );
        }
    }
}