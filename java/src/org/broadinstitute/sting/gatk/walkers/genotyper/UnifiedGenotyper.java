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
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotator;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.cmdLine.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.geli.GeliGenotypeWriter;
import org.broadinstitute.sting.utils.genotype.glf.GLFGenotypeWriter;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.util.*;
import java.io.PrintStream;


/**
 * A variant caller which unifies the approaches of several disparate callers.  Works for single-sample,
 * multi-sample, and pooled data.  The user can choose from several different incorporated calculation models.
 */
@Reference(window=@Window(start=-20,stop=20))
public class UnifiedGenotyper extends LocusWalker<VariantCallContext, UnifiedGenotyper.UGStatistics> implements TreeReducible<UnifiedGenotyper.UGStatistics> {
    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Argument(doc = "File to which variants should be written", required = false)
    public GenotypeWriter writer = null;

    @Argument(fullName = "verbose_mode", shortName = "verbose", doc = "File to print all of the annotated and detailed debugging output", required = false)
    public PrintStream verboseWriter = null;

    @Argument(fullName = "beagle_file", shortName = "beagle", doc = "File to print BEAGLE-specific data for use with imputation", required = false)
    public PrintStream beagleWriter = null;

    // the model used for calculating genotypes
    private ThreadLocal<GenotypeCalculationModel> gcm = new ThreadLocal<GenotypeCalculationModel>();

    // samples in input
    private Set<String> samples = new HashSet<String>();

    // should we annotate dbsnp?
    private boolean annotateDbsnp = false;

    // Enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /**
     * Inner class for collecting output statistics from the UG
     */
    public class UGStatistics {
        /** The total number of passes examined -- i.e., the number of map calls */
        long nBasesVisited = 0;

        /** The number of bases that were potentially callable -- i.e., those not at excessive coverage or masked with N */
        long nBasesCallable = 0;

        /** The number of bases called confidently (according to user threshold), either ref or other */
        long nBasesCalledConfidently = 0;

        double percentCallableOfAll()    { return (100.0 * nBasesCallable) / nBasesVisited; }
        double percentCalledOfAll()      { return (100.0 * nBasesCalledConfidently) / nBasesVisited; }
        double percentCalledOfCallable() { return (100.0 * nBasesCalledConfidently) / nBasesCallable; }
    }

    /**
     * Sets the argument collection for the UnifiedGenotyper.
     * To be used with walkers that call the UnifiedGenotyper's map function
     * and consequently can't set these arguments on the command-line
     *
     * @param UAC the UnifiedArgumentCollection
     *
     **/
    public void setUnifiedArgumentCollection(UnifiedArgumentCollection UAC) {
        this.UAC = UAC;
        initialize();
    }

    /**
     * Initialize the samples, output, and genotype calculation model
     *
     **/
    public void initialize() {
        // deal with input errors
        if ( UAC.POOLSIZE > 0 && UAC.genotypeModel != GenotypeCalculationModel.Model.POOLED ) {
            throw new IllegalArgumentException("Attempting to use a model other than POOLED with pooled data. Please set the model to POOLED.");
        }
        if ( UAC.POOLSIZE < 1 && UAC.genotypeModel == GenotypeCalculationModel.Model.POOLED ) {
            throw new IllegalArgumentException("Attempting to use the POOLED model with a pool size less than 1. Please set the pool size to an appropriate value.");
        }
        if ( beagleWriter != null && UAC.genotypeModel == GenotypeCalculationModel.Model.EM_POINT_ESTIMATE ) {
            throw new IllegalArgumentException("BEAGLE output is not currently supported in the EM_POINT_ESTIMATE calculation model.");
        }
        if ( getToolkit().getArguments().numberOfThreads > 1 && UAC.ASSUME_SINGLE_SAMPLE != null ) {
            // the ASSUME_SINGLE_SAMPLE argument can't be handled (at least for now) while we are multi-threaded because the IO system doesn't know how to get the sample name
            throw new IllegalArgumentException("For technical reasons, the ASSUME_SINGLE_SAMPLE argument cannot be used with multiple threads");
        }

        // get all of the unique sample names - unless we're in POOLED mode, in which case we ignore the sample names
        if ( UAC.genotypeModel != GenotypeCalculationModel.Model.POOLED ) {
            // if we're supposed to assume a single sample, do so
            if ( UAC.ASSUME_SINGLE_SAMPLE != null )
                samples.add(UAC.ASSUME_SINGLE_SAMPLE);
            else
                samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());

            // for ( String sample : samples )
            //     logger.debug("SAMPLE: " + sample);
        }

        // in pooled mode we need to check that the format is acceptable
        if ( UAC.genotypeModel == GenotypeCalculationModel.Model.POOLED && writer != null ) {
            // only multi-sample calls use Variations
            if ( !writer.supportsMultiSample() )
                throw new IllegalArgumentException("The POOLED model is not compatible with the specified format; try using VCF instead");

            // when using VCF with multiple threads, we need to turn down the validation stringency so that writing temporary files will work
            if ( getToolkit().getArguments().numberOfThreads > 1 && writer instanceof VCFGenotypeWriter )
                ((VCFGenotypeWriter)writer).setValidationStringency(VCFGenotypeWriterAdapter.VALIDATION_STRINGENCY.SILENT);
        }

        // initialize the writers
        if ( verboseWriter != null ) {
            if(UAC.genotypeModel != GenotypeCalculationModel.Model.EM_POINT_ESTIMATE) {
                StringBuilder header = new StringBuilder("AFINFO\tLOC\tMAF\tF\tNullAFpriors\t");
                for ( char altAllele : BaseUtils.BASES ) {
                    char base = Character.toUpperCase(altAllele);
                    header.append("POfDGivenAFFor" + base + "\t");
                    header.append("PosteriorAFFor" + base + "\t");
                }
                verboseWriter.println(header);
            }
        }
        if ( beagleWriter != null ) {
            beagleWriter.print("marker alleleA alleleB");
            for ( String sample : samples )
                beagleWriter.print(String.format(" %s %s %s", sample, sample, sample));
            beagleWriter.println();
        }

        // check to see whether a dbsnp rod was included
        List<ReferenceOrderedDataSource> dataSources = getToolkit().getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            ReferenceOrderedData rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(rodDbSNP.class) ) {
                annotateDbsnp = true;
                break;
            }
        }

        // *** If we were called by another walker, then we don't ***
        // *** want to do any of the other initialization steps.  ***
        if ( writer == null )
            return;

        // *** If we got here, then we were instantiated by the GATK engine ***

        // initialize the header
        GenotypeWriterFactory.writeHeader(writer, GenomeAnalysisEngine.instance.getSAMFileHeader(), samples, getHeaderInfo());
    }

    private Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // this is only applicable to VCF
        if ( !(writer instanceof VCFGenotypeWriter) )
            return headerInfo;

        // first, the basic info
        headerInfo.add(new VCFHeaderLine("source", "UnifiedGenotyper"));
        headerInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        // annotation (INFO) fields from VariantAnnotator
        if ( UAC.ALL_ANNOTATIONS )
            headerInfo.addAll(VariantAnnotator.getAllVCFAnnotationDescriptions());
        else
            headerInfo.addAll(VariantAnnotator.getVCFAnnotationDescriptions());

        // annotation (INFO) fields from UnifiedGenotyper
        headerInfo.add(new VCFInfoHeaderLine(VCFRecord.ALLELE_FREQUENCY_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Allele Frequency"));
        if ( annotateDbsnp )
            headerInfo.add(new VCFInfoHeaderLine(VCFRecord.DBSNP_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "dbSNP membership"));
        if ( !UAC.NO_SLOD )
            headerInfo.add(new VCFInfoHeaderLine(VCFRecord.STRAND_BIAS_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Strand Bias"));

        // FORMAT fields if not in POOLED mode
        if ( UAC.genotypeModel != GenotypeCalculationModel.Model.POOLED ) {
            headerInfo.addAll(VCFGenotypeRecord.getSupportedHeaderStrings());
            headerInfo.add(new VCFInfoHeaderLine(VCFRecord.ALLELE_COUNT_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"));
            headerInfo.add(new VCFInfoHeaderLine(VCFRecord.ALLELE_NUMBER_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Total number of alleles in called genotypes"));
        }

        // all of the arguments from the argument collection
        Set<Object> args = new HashSet<Object>();
        args.add(UAC);
        args.addAll(getToolkit().getFilters());
        Map<String,String> commandLineArgs = CommandLineUtils.getApproximateCommandLineArguments(args);
        for ( Map.Entry<String, String> commandLineArg : commandLineArgs.entrySet() )
            headerInfo.add(new VCFHeaderLine(String.format("UG_%s", commandLineArg.getKey()), commandLineArg.getValue()));            

        return headerInfo;
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     */
    public VariantCallContext map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {

        // initialize the GenotypeCalculationModel for this thread if that hasn't been done yet
        if ( gcm.get() == null ) {
            GenotypeWriterFactory.GENOTYPE_FORMAT format = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;
            if ( writer != null ) {
                if ( writer instanceof VCFGenotypeWriter )
                    format = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;
                else if ( writer instanceof GLFGenotypeWriter )
                    format = GenotypeWriterFactory.GENOTYPE_FORMAT.GLF;
                else if ( writer instanceof GeliGenotypeWriter )
                    format = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;
                else
                    throw new StingException("Unsupported genotype format: " + writer.getClass().getName());
            }
            gcm.set(GenotypeCalculationModelFactory.makeGenotypeCalculation(samples, logger, UAC, format, verboseWriter, beagleWriter));
        }

        char ref = Character.toUpperCase(refContext.getBase());
        if ( !BaseUtils.isRegularBase(ref) )
            return null;

        // don't try to call if we couldn't read in all reads at this locus (since it wasn't properly downsampled)
        if ( rawContext.hasExceededMaxPileup() )
            return null;

        ReadBackedPileup rawPileup = rawContext.getBasePileup();
        
        // filter the context based on min base and mapping qualities
        ReadBackedPileup pileup = rawPileup.getBaseAndMappingFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE, UAC.MIN_MAPPING_QUALTY_SCORE);

        // filter the context based on mapping quality and mismatch rate
        pileup = filterPileup(pileup, refContext, UAC);

        // don't call when there is no coverage
        if ( pileup.size() == 0 )
            return null;

        // are there too many deletions in the pileup?
        if ( isValidDeletionFraction(UAC.MAX_DELETION_FRACTION) &&
             (double)pileup.getNumberOfDeletions() / (double)pileup.size() > UAC.MAX_DELETION_FRACTION )
            return null;

        // stratify the AlignmentContext and cut by sample
        // Note that for testing purposes, we may want to throw multi-samples at pooled mode
        Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(pileup, UAC.ASSUME_SINGLE_SAMPLE, (UAC.genotypeModel == GenotypeCalculationModel.Model.POOLED ? PooledCalculationModel.POOL_SAMPLE_NAME : null));
        if ( stratifiedContexts == null )
            return null;

        DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, UAC.heterozygosity, DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);
        VariantCallContext call = gcm.get().callLocus(tracker, ref, rawContext.getLocation(), stratifiedContexts, priors);

        // annotate the call, if possible
        if ( call != null && call.variation != null && call.variation instanceof ArbitraryFieldsBacked ) {
            // first off, we want to use the *unfiltered* context for the annotations
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(rawContext.getBasePileup());

            Map<String, String> annotations;
            if ( UAC.ALL_ANNOTATIONS )
                annotations = VariantAnnotator.getAllAnnotations(tracker, refContext, stratifiedContexts, call.variation, annotateDbsnp);
            else
                annotations = VariantAnnotator.getAnnotations(tracker, refContext, stratifiedContexts, call.variation, annotateDbsnp);
            ((ArbitraryFieldsBacked)call.variation).setFields(annotations);
        }

        return call;
    }

    // filter based on maximum mismatches and bad mates
    private static ReadBackedPileup filterPileup(ReadBackedPileup pileup, ReferenceContext refContext, UnifiedArgumentCollection UAC) {

        ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();
        for ( PileupElement p : pileup ) {
            if  ( (UAC.USE_BADLY_MATED_READS || !p.getRead().getReadPairedFlag() || p.getRead().getMateUnmappedFlag() || p.getRead().getMateReferenceIndex() == p.getRead().getReferenceIndex()) &&
                  AlignmentUtils.mismatchesInRefWindow(p, refContext, true) <= UAC.MAX_MISMATCHES )
                filteredPileup.add(p);
        }
        return new ReadBackedPileup(pileup.getLocation(), filteredPileup);
    }

    private static boolean isValidDeletionFraction(double d) {
        return ( d >= 0.0 && d <= 1.0 );
    }

    // ------------------------------------------------------------------------------------------------
    //
    // Reduce
    //
    // ------------------------------------------------------------------------------------------------
    public UGStatistics reduceInit() { return new UGStatistics(); }

    public UGStatistics treeReduce(UGStatistics lhs, UGStatistics rhs) {
        lhs.nBasesCallable += rhs.nBasesCallable;
        lhs.nBasesCalledConfidently += rhs.nBasesCalledConfidently;
        lhs.nBasesVisited += rhs.nBasesVisited;
        return lhs;
    }

    public UGStatistics reduce(VariantCallContext value, UGStatistics sum) {
        // We get a point for reaching reduce :-)
        sum.nBasesVisited++;

        // can't call the locus because of no coverage
        if ( value == null )
            return sum;

        // A call was attempted -- the base was potentially callable
        sum.nBasesCallable++;

        // if the base was confidently called something, print it out 
        sum.nBasesCalledConfidently += value.confidentlyCalled ? 1 : 0;

        // can't make a confident variant call here
        if ( value.genotypes == null ||
                (UAC.genotypeModel != GenotypeCalculationModel.Model.POOLED && value.genotypes.size() == 0) ) {
            return sum;
        }

        // if we have a single-sample call (single sample from PointEstimate model returns no VariationCall data)
        if ( value.variation == null || (!writer.supportsMultiSample() && samples.size() <= 1) ) {
            writer.addGenotypeCall(value.genotypes.get(0));
        }

        // use multi-sample mode if we have multiple samples or the output type allows it
        else {
            writer.addMultiSampleCall(value.genotypes, value.variation);
        }

        return sum;
    }

    // Close any file writers
    public void onTraversalDone(UGStatistics sum) {
        logger.info(String.format("Visited bases                                %d", sum.nBasesVisited));
        logger.info(String.format("Callable bases                               %d", sum.nBasesCallable));
        logger.info(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
        logger.info(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
        logger.info(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
        logger.info(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
//        logger.info("Processed " + sum.nBasesCallable + " loci that are callable for SNPs");
    }
}
