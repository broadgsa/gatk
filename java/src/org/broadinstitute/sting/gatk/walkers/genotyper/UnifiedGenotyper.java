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

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotator;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.ArgumentCollection;
import org.broadinstitute.sting.utils.genotype.*;

import java.io.File;
import java.util.*;


@Reference(window=@Window(start=-20,stop=20))
public class UnifiedGenotyper extends LocusWalker<Pair<List<Genotype>, GenotypeLocusData>, Integer> {

    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = false)
    public File VARIANTS_FILE = null;


    // the model used for calculating genotypes
    private GenotypeCalculationModel gcm;

    // output writer
    private GenotypeWriter writer;

    // samples in input
    private HashSet<String> samples;

    // keep track of some metrics about our calls
    private CallMetrics callsMetrics;


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
            throw new StingException("Attempting to use a model other than POOLED with pooled data. Please set the model to POOLED.");
        }
        if ( UAC.LOD_THRESHOLD > Double.MIN_VALUE ) {
            StringBuilder sb = new StringBuilder();
            sb.append("\n***\tThe --lod_threshold argument is no longer supported; instead, please use --min_confidence_threshold.");
            sb.append("\n***\tThere is approximately a 10-to-1 mapping from confidence to LOD.");
            sb.append("\n***\tUse Q" + (10.0 * UAC.LOD_THRESHOLD) + " as an approximate equivalent to your LOD " + UAC.LOD_THRESHOLD + " cutoff");
            throw new StingException(sb.toString());
        }

        // get all of the unique sample names
        samples = new HashSet<String>();
        // if we're supposed to assume a single sample
        if ( UAC.ASSUME_SINGLE_SAMPLE != null ) {
            samples.add(UAC.ASSUME_SINGLE_SAMPLE);
        } else {
            List<SAMReadGroupRecord> readGroups = getToolkit().getSAMFileHeader().getReadGroups();
            for ( SAMReadGroupRecord readGroup : readGroups )
                samples.add(readGroup.getSample());
        }

        // print them out for debugging (need separate loop to ensure uniqueness)
        for ( String sample : samples )
            logger.debug("SAMPLE: " + sample);

        gcm = GenotypeCalculationModelFactory.makeGenotypeCalculation(samples, logger, UAC);

        // *** If we were called by another walker, then we don't ***
        // *** want to do any of the other initialization steps.  ***
        if ( VARIANTS_FILE == null && out == null )
            return;

        // if we got here, then we were instantiated by the GATK engine

        // create the output writer stream
        if ( VARIANTS_FILE != null )
            writer = GenotypeWriterFactory.create(UAC.VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), VARIANTS_FILE,
                                                  "UnifiedGenotyper",
                                                  this.getToolkit().getArguments().referenceFile.getName(),
                                                  samples);
        else
            writer = GenotypeWriterFactory.create(UAC.VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), out,
                                                  "UnifiedGenotyper",
                                                  this.getToolkit().getArguments().referenceFile.getName(),
                                                  samples);
        callsMetrics = new CallMetrics();
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param fullContext contextual information around the locus
     */
    public Pair<List<Genotype>, GenotypeLocusData> map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext fullContext) {
        char ref = Character.toUpperCase(refContext.getBase());
        if ( !BaseUtils.isRegularBase(ref) )
            return null;

        // remove mapping quality zero reads
        AlignmentContext MQ0freeContext = filterAlignmentContext(fullContext);

        // an optimization to speed things up when there is no coverage or when overly covered
        if ( MQ0freeContext.getReads().size() == 0 ||
             (UAC.MAX_READS_IN_PILEUP > 0 && MQ0freeContext.getReads().size() > UAC.MAX_READS_IN_PILEUP) )
            return null;

        DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, UAC.heterozygosity, DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);
        Pair<List<Genotype>, GenotypeLocusData> call = gcm.calculateGenotype(tracker, ref, MQ0freeContext, priors);

        // annotate the call, if possible
        if ( call != null && call.second != null && call.second instanceof ArbitraryFieldsBacked && ! (VARIANTS_FILE == null && out == null) ) {
            Map<String, String> annotations = VariantAnnotator.getAnnotations(refContext, fullContext, call.first);
            ((ArbitraryFieldsBacked)call.second).setFields(annotations);
        }

        return call;
    }

    private AlignmentContext filterAlignmentContext(AlignmentContext context) {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        List<SAMRecord> newReads = new ArrayList<SAMRecord>();
        List<Integer> newOffsets = new ArrayList<Integer>();

        for (int i = 0; i < reads.size(); i++) {
            SAMRecord read = reads.get(i);
            if ( read.getMappingQuality() != 0 ) {
                newReads.add(read);
                newOffsets.add(offsets.get(i));
            }
        }

        return new AlignmentContext(context.getLocation(), newReads, newOffsets);                
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Pair<List<Genotype>, GenotypeLocusData> value, Integer sum) {
        // can't call the locus because of no coverage
        if ( value == null )
            return sum;

        callsMetrics.nCalledBases++;

        // can't make a confident variant call here
        if ( value.first == null || value.first.size() == 0 ) {
            callsMetrics.nNonConfidentCalls++;
            return sum;
        }

        callsMetrics.nConfidentCalls++;

        // if we have a single-sample call (single sample from PointEstimate model returns no genotype locus data)
        if ( value.second == null || (!writer.supportsMultiSample() && samples.size() == 1) ) {
            writer.addGenotypeCall(value.first.get(0));
        }

        // use multi-sample mode if we have multiple samples or the output type allows it
        else {
            writer.addMultiSampleCall(value.first, value.second);
        }

        return sum + 1;
    }

    // Close any file writers
    public void onTraversalDone(Integer sum) {
        writer.close();
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
