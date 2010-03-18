/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.MutableVariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.geli.GeliGenotypeWriter;
import org.broadinstitute.sting.utils.genotype.glf.GLFGenotypeWriter;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.PrintStream;


public class UnifiedGenotyperEngine {

    // should we annotate dbsnp?
    protected boolean annotateDbsnp = false;
    // should we annotate hapmap2?
    protected boolean annotateHapmap2 = false;
    // should we annotate hapmap3?
    protected boolean annotateHapmap3 = false;

    // the unified argument collection
    protected UnifiedArgumentCollection UAC = null;

    // the annotation engine
    protected VariantAnnotatorEngine annotationEngine;

    // the model used for calculating genotypes
    protected ThreadLocal<GenotypeCalculationModel> gcm = new ThreadLocal<GenotypeCalculationModel>();

    // the various loggers and writers
    protected Logger logger = null;
    protected GenotypeWriter genotypeWriter = null;
    protected PrintStream verboseWriter = null;
    protected PrintStream beagleWriter = null;

    // samples in input
    protected Set<String> samples = new HashSet<String>();


    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC) {
        initialize(toolkit, UAC, null, null, null, null, null);
    }

    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, GenotypeWriter genotypeWriter, PrintStream verboseWriter, PrintStream beagleWriter, VariantAnnotatorEngine engine) {
        initialize(toolkit, UAC, logger, genotypeWriter, verboseWriter, beagleWriter, engine);

    }

    private void initialize(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, GenotypeWriter genotypeWriter, PrintStream verboseWriter, PrintStream beagleWriter, VariantAnnotatorEngine engine) {
        this.UAC = UAC;
        this.logger = logger;
        this.genotypeWriter = genotypeWriter;
        this.verboseWriter = verboseWriter;
        this.beagleWriter = beagleWriter;
        this.annotationEngine = engine;

        // deal with input errors
        if ( UAC.genotypeModel == GenotypeCalculationModel.Model.INDELS && !(genotypeWriter instanceof VCFGenotypeWriter) ) {
            throw new IllegalArgumentException("Attempting to use an output format other than VCF with indels. Please set the output format to VCF.");
        }
        if ( UAC.POOLSIZE > 0 && UAC.genotypeModel != GenotypeCalculationModel.Model.POOLED ) {
            throw new IllegalArgumentException("Attempting to use a model other than POOLED with pooled data. Please set the model to POOLED.");
        }
        if ( UAC.POOLSIZE < 1 && UAC.genotypeModel == GenotypeCalculationModel.Model.POOLED ) {
            throw new IllegalArgumentException("Attempting to use the POOLED model with a pool size less than 1. Please set the pool size to an appropriate value.");
        }
        if ( toolkit.getArguments().numberOfThreads > 1 && UAC.ASSUME_SINGLE_SAMPLE != null ) {
            // the ASSUME_SINGLE_SAMPLE argument can't be handled (at least for now) while we are multi-threaded because the IO system doesn't know how to get the sample name
            throw new IllegalArgumentException("For technical reasons, the ASSUME_SINGLE_SAMPLE argument cannot be used with multiple threads");
        }

        // get all of the unique sample names - unless we're in POOLED mode, in which case we ignore the sample names
        if ( UAC.genotypeModel != GenotypeCalculationModel.Model.POOLED ) {
            // if we're supposed to assume a single sample, do so
            if ( UAC.ASSUME_SINGLE_SAMPLE != null )
                this.samples.add(UAC.ASSUME_SINGLE_SAMPLE);
            else
                this.samples = SampleUtils.getSAMFileSamples(toolkit.getSAMFileHeader());
        }

        // in pooled mode we need to check that the format is acceptable
        if ( UAC.genotypeModel == GenotypeCalculationModel.Model.POOLED && genotypeWriter != null ) {
            // when using VCF with multiple threads, we need to turn down the validation stringency so that writing temporary files will work
            if ( toolkit.getArguments().numberOfThreads > 1 && genotypeWriter instanceof VCFGenotypeWriter )
                ((VCFGenotypeWriter)genotypeWriter).setValidationStringency(VCFGenotypeWriterAdapter.VALIDATION_STRINGENCY.SILENT);
        }

        // check to see whether a dbsnp rod was included
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            ReferenceOrderedData rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(rodDbSNP.class) ) {
                this.annotateDbsnp = true;
            }
            if ( rod.getName().equals("hapmap2") ) {
                this.annotateHapmap2 = true;
            }
            if ( rod.getName().equals("hapmap3") ) {
                this.annotateHapmap3 = true;
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
            GenotypeWriterFactory.GENOTYPE_FORMAT format = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;
            if ( genotypeWriter != null ) {
                if ( genotypeWriter instanceof VCFGenotypeWriter )
                    format = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;
                else if ( genotypeWriter instanceof GLFGenotypeWriter)
                    format = GenotypeWriterFactory.GENOTYPE_FORMAT.GLF;
                else if ( genotypeWriter instanceof GeliGenotypeWriter)
                    format = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;
                else
                    throw new StingException("Unsupported genotype format: " + genotypeWriter.getClass().getName());
            }
            gcm.set(GenotypeCalculationModelFactory.makeGenotypeCalculation(samples, logger, UAC, format, verboseWriter, beagleWriter));
        }

        char ref = Character.toUpperCase(refContext.getBase());
        if ( !BaseUtils.isRegularBase(ref) )
            return null;

        // don't try to call if we couldn't read in all reads at this locus (since it wasn't properly downsampled)
        if ( rawContext.hasExceededMaxPileup() )
            return null;

        VariantCallContext call;

        if ( rawContext.hasExtendedEventPileup() ) {

            ReadBackedExtendedEventPileup rawPileup = rawContext.getExtendedEventPileup();

            // filter the context based on min mapping quality
            ReadBackedExtendedEventPileup pileup = rawPileup.getMappingFilteredPileup(UAC.MIN_MAPPING_QUALTY_SCORE);

            // filter the context based on bad mates and mismatch rate
            pileup = filterPileup(pileup, refContext);

            // don't call when there is no coverage
            if ( pileup.size() == 0 )
                return null;

            // stratify the AlignmentContext and cut by sample
            Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(pileup, UAC.ASSUME_SINGLE_SAMPLE, null);
            if ( stratifiedContexts == null )
                return null;

            call = gcm.get().callExtendedLocus(tracker, ref, rawContext.getLocation(), stratifiedContexts);

        } else {

            ReadBackedPileup rawPileup = rawContext.getBasePileup();

            // filter the context based on min base and mapping qualities
            ReadBackedPileup pileup = rawPileup.getBaseAndMappingFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE, UAC.MIN_MAPPING_QUALTY_SCORE);

            // filter the context based on bad mates and mismatch rate
            pileup = filterPileup(pileup, refContext);

            // don't call when there is no coverage
            if ( pileup.size() == 0 )
                return null;

            // are there too many deletions in the pileup?
            if ( UAC.genotypeModel != GenotypeCalculationModel.Model.INDELS &&
                 isValidDeletionFraction(UAC.MAX_DELETION_FRACTION) &&
                 (double)pileup.getNumberOfDeletions() / (double)pileup.size() > UAC.MAX_DELETION_FRACTION )
                return null;

            // stratify the AlignmentContext and cut by sample
            // Note that for testing purposes, we may want to throw multi-samples at pooled mode
            Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(pileup, UAC.ASSUME_SINGLE_SAMPLE, (UAC.genotypeModel == GenotypeCalculationModel.Model.POOLED ? PooledCalculationModel.POOL_SAMPLE_NAME : null));
            if ( stratifiedContexts == null )
                return null;

            DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, UAC.heterozygosity, DiploidGenotypePriors.PROB_OF_REFERENCE_ERROR);
            call = gcm.get().callLocus(tracker, ref, rawContext.getLocation(), stratifiedContexts, priors);

            // annotate the call, if possible
            if ( call != null && call.vc != null ) {
                // first off, we want to use the *unfiltered* context for the annotations
                stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(rawContext.getBasePileup());
                annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, (MutableVariantContext)call.vc);
            }
        }

        return call;
    }

    // filter based on maximum mismatches and bad mates
    private ReadBackedPileup filterPileup(ReadBackedPileup pileup, ReferenceContext refContext) {

        ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();
        for ( PileupElement p : pileup ) {
            if  ( (UAC.USE_BADLY_MATED_READS || !p.getRead().getReadPairedFlag() || p.getRead().getMateUnmappedFlag() || p.getRead().getMateReferenceIndex() == p.getRead().getReferenceIndex()) &&
                  AlignmentUtils.mismatchesInRefWindow(p, refContext, true) <= UAC.MAX_MISMATCHES )
                filteredPileup.add(p);
        }
        return new ReadBackedPileup(pileup.getLocation(), filteredPileup);
    }

    // filter based on maximum mismatches and bad mates
    private ReadBackedExtendedEventPileup filterPileup(ReadBackedExtendedEventPileup pileup, ReferenceContext refContext) {

        ArrayList<ExtendedEventPileupElement> filteredPileup = new ArrayList<ExtendedEventPileupElement>();
        for ( ExtendedEventPileupElement p : pileup ) {
            if  ( (UAC.USE_BADLY_MATED_READS || !p.getRead().getReadPairedFlag() || p.getRead().getMateUnmappedFlag() || p.getRead().getMateReferenceIndex() == p.getRead().getReferenceIndex()) &&
                  AlignmentUtils.mismatchesInRefWindow(p, refContext, true) <= UAC.MAX_MISMATCHES )
                filteredPileup.add(p);
        }
        return new ReadBackedExtendedEventPileup(pileup.getLocation(), filteredPileup);
    }

    private static boolean isValidDeletionFraction(double d) {
        return ( d >= 0.0 && d <= 1.0 );
    }
}