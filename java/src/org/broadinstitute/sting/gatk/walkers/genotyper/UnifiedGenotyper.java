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
import org.broadinstitute.sting.gatk.filters.*;
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


/**
 * A variant caller which unifies the approaches of several disparate callers.  Works for single-sample,
 * multi-sample, and pooled data.  The user can choose from several different incorporated calculation models.
 */
@Reference(window=@Window(start=-20,stop=20))
@ReadFilters({ZeroMappingQualityReadFilter.class,MappingQualityReadFilter.class,BadMateReadFilter.class})
public class UnifiedGenotyper extends LocusWalker<Pair<VariationCall, List<Genotype>>, Integer> implements TreeReducible<Integer> {

    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Argument(doc = "File to which variants should be written", required = false)
    public GenotypeWriter writer = null;

    // the model used for calculating genotypes
    private GenotypeCalculationModel gcm;

    // samples in input
    private Set<String> samples;


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
        gcm.close();
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
        if ( UAC.LOD_THRESHOLD > Double.MIN_VALUE ) {
            StringBuilder sb = new StringBuilder();
            sb.append("\n***\tThe --lod_threshold argument is no longer supported; instead, please use --min_confidence_threshold.");
            sb.append("\n***\tThere is approximately a 10-to-1 mapping from confidence to LOD.");
            sb.append("\n***\tUse Q" + (10.0 * UAC.LOD_THRESHOLD) + " as an approximate equivalent to your LOD " + UAC.LOD_THRESHOLD + " cutoff");
            throw new IllegalArgumentException(sb.toString());
        }

        // get all of the unique sample names
        // if we're supposed to assume a single sample
        if ( UAC.ASSUME_SINGLE_SAMPLE != null )
            samples.add(UAC.ASSUME_SINGLE_SAMPLE);
        else
            samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());

        // print them out for debugging (need separate loop to ensure uniqueness)
        // for ( String sample : samples )
        //     logger.debug("SAMPLE: " + sample);

        GenotypeWriterFactory.GENOTYPE_FORMAT format = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;
        if(writer != null) {
            if(writer instanceof VCFGenotypeWriter)
                format = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;
            else if(writer instanceof GLFGenotypeWriter)
                format = GenotypeWriterFactory.GENOTYPE_FORMAT.GLF;
            else if(writer instanceof GeliGenotypeWriter)
                format = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;
            else
                throw new StingException("Unsupported genotype format: " + writer.getClass().getName());
        }

        gcm = GenotypeCalculationModelFactory.makeGenotypeCalculation(samples, logger, UAC, format);

        // *** If we were called by another walker, then we don't ***
        // *** want to do any of the other initialization steps.  ***
        if ( writer == null )
            return;

        // *** If we got here, then we were instantiated by the GATK engine ***

        // if we're in POOLED mode, we no longer want the sample names
        // (although they're still stored in the gcm if needed)
        if ( UAC.genotypeModel == GenotypeCalculationModel.Model.POOLED )
            samples.clear();

        // get the optional header fields
        Set<VCFHeaderLine> headerInfo = getHeaderInfo();

        // initialize the header
        GenotypeWriterFactory.writeHeader(writer, GenomeAnalysisEngine.instance.getSAMFileHeader(), samples, headerInfo);
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
        headerInfo.add(new VCFHeaderLine("INFO_NOTE", "\"All annotations in the INFO field are generated only from the FILTERED context used for calling variants\""));
        headerInfo.add(new VCFInfoHeaderLine("AF", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Allele Frequency"));
        headerInfo.add(new VCFInfoHeaderLine("NS", 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Number of Samples With Data"));
        if ( !UAC.NO_SLOD )
            headerInfo.add(new VCFInfoHeaderLine("SB", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Strand Bias"));

        // FORMAT fields if not in POOLED mode
        if ( UAC.genotypeModel != GenotypeCalculationModel.Model.POOLED )
            headerInfo.addAll(VCFGenotypeRecord.getSupportedHeaderStrings());

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
        char ref = Character.toUpperCase(refContext.getBase());
        if ( !BaseUtils.isRegularBase(ref) )
            return null;

        // filter the context based on min base and mapping qualities
        ReadBackedPileup pileup = rawContext.getPileup().getBaseFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE);

        // filter the context based on mismatches
        pileup = filterPileup(pileup, refContext, UAC.MAX_MISMATCHES);

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
        Pair<VariationCall, List<Genotype>> call = gcm.calculateGenotype(tracker, ref, rawContext.getLocation(), stratifiedContexts, priors);

        // annotate the call, if possible
        if ( call != null && call.first != null && call.first instanceof ArbitraryFieldsBacked ) {
            Map<String, String> annotations;
            if ( UAC.ALL_ANNOTATIONS )
                annotations = VariantAnnotator.getAllAnnotations(refContext, stratifiedContexts, call.first);
            else
                annotations = VariantAnnotator.getAnnotations(refContext, stratifiedContexts, call.first);
            ((ArbitraryFieldsBacked)call.first).setFields(annotations);
        }

        return call;
    }

    // filter based on maximum mismatches and bad mates
    private static ReadBackedPileup filterPileup(ReadBackedPileup pileup, ReferenceContext refContext, int maxMismatches) {
        ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();
        for ( PileupElement p : pileup ) {
            if ( AlignmentUtils.mismatchesInRefWindow(p, refContext, true) <= maxMismatches )
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
        gcm.close();
        logger.info("Processed " + sum + " loci that are callable for SNPs");
    }

    /**
     * A class to keep track of some basic metrics about our calls
     */
    protected class CallMetrics {
        long nConfidentCalls = 0;
        long nNonConfidentCalls = 0;
        long nCalledBases = 0;

        CallMetrics() {}

        public String toString() {
            return String.format("UG: %d confident and %d non-confident calls were made at %d bases",
                    nConfidentCalls, nNonConfidentCalls, nCalledBases);
        }
    }      
}
