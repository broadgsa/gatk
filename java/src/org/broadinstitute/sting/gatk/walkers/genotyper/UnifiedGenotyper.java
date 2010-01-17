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
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
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
public class UnifiedGenotyper extends LocusWalker<Pair<VariationCall, List<Genotype>>, Integer> implements TreeReducible<Integer> {

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


    /** Enable deletions in the pileup **/
    public boolean includeReadsWithDeletionAtLoci() { return true; }

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
                    char base = Character.toLowerCase(altAllele);
                    header.append("POfDGivenAFFor" + base + "\t");
                    header.append("PosteriorAFFor" + base + "\t");
                }
                verboseWriter.println(header);
            }
        }
        if ( beagleWriter != null ) {
            beagleWriter.print("marker alleleA alleleB");
            for ( String sample : samples ) {
                beagleWriter.print(' ');
                beagleWriter.print(sample);
            }
            beagleWriter.println();
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
    public Pair<VariationCall, List<Genotype>> map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {

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

        // filter the context based on min base and mapping qualities
        ReadBackedPileup pileup = rawContext.getBasePileup().getBaseFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE);

        // filter the context based on mapping quality and mismatch rate
        pileup = filterPileup(pileup, refContext, UAC);

        // an optimization to speed things up when there is no coverage or when overly covered
        if ( pileup.size() == 0 ||
             (UAC.MAX_READS_IN_PILEUP > 0 && pileup.size() > UAC.MAX_READS_IN_PILEUP) )
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
        Pair<VariationCall, List<Genotype>> call = gcm.get().calculateGenotype(tracker, ref, rawContext.getLocation(), stratifiedContexts, priors);

        // annotate the call, if possible
        if ( call != null && call.first != null && call.first instanceof ArbitraryFieldsBacked ) {
            // first off, we want to use the *unfiltered* context for the annotations
            stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(rawContext.getBasePileup());

            Map<String, String> annotations;
            if ( UAC.ALL_ANNOTATIONS )
                annotations = VariantAnnotator.getAllAnnotations(tracker, refContext, stratifiedContexts, call.first);
            else
                annotations = VariantAnnotator.getAnnotations(tracker, refContext, stratifiedContexts, call.first);
            ((ArbitraryFieldsBacked)call.first).setFields(annotations);
        }

        return call;
    }

    // filter based on maximum mismatches and bad mates
    private static ReadBackedPileup filterPileup(ReadBackedPileup pileup, ReferenceContext refContext, UnifiedArgumentCollection UAC) {

        ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();
        for ( PileupElement p : pileup ) {
            if  ( p.getMappingQual() >= UAC.MIN_MAPPING_QUALTY_SCORE &&
                  (UAC.USE_BADLY_MATED_READS || !p.getRead().getReadPairedFlag() || p.getRead().getMateUnmappedFlag() || p.getRead().getMateReferenceIndex() == p.getRead().getReferenceIndex()) &&
                  AlignmentUtils.mismatchesInRefWindow(p, refContext, true) <= UAC.MAX_MISMATCHES )
                filteredPileup.add(p);
        }
        return new ReadBackedPileup(pileup.getLocation(), filteredPileup);
    }

    private static boolean isValidDeletionFraction(double d) {
        return ( d >= 0.0 && d <= 1.0 );
    }

    public Integer reduceInit() { return 0; }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;        
    }

    public Integer reduce(Pair<VariationCall, List<Genotype>> value, Integer sum) {
        // can't call the locus because of no coverage
        if ( value == null )
            return sum;

        // can't make a confident variant call here
        if ( value.second == null ||
                (UAC.genotypeModel != GenotypeCalculationModel.Model.POOLED && value.second.size() == 0) ) {
            return sum;
        }

        // if we have a single-sample call (single sample from PointEstimate model returns no VariationCall data)
        if ( value.first == null || (!writer.supportsMultiSample() && samples.size() <= 1) ) {
            writer.addGenotypeCall(value.second.get(0));
        }

        // use multi-sample mode if we have multiple samples or the output type allows it
        else {
            writer.addMultiSampleCall(value.second, value.first);
        }

        return sum + 1;
    }

    // Close any file writers
    public void onTraversalDone(Integer sum) {
        logger.info("Processed " + sum + " loci that are callable for SNPs");
    }
}
