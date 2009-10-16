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
import org.broadinstitute.sting.gatk.filters.MissingReadGroupFilter;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.ArgumentCollection;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.GenotypeMetaData;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;


@ReadFilters({ZeroMappingQualityReadFilter.class, MissingReadGroupFilter.class})
public class UnifiedGenotyper extends LocusWalker<Pair<List<GenotypeCall>, GenotypeMetaData>, Integer> {

    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = false)
    public File VARIANTS_FILE = null;

    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File format to be used; default is VCF", required = false)
    public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;


    // the model used for calculating genotypes
    private GenotypeCalculationModel gcm;

    // output writer
    private GenotypeWriter writer;

    // samples in input
    private HashSet<String> samples;

    // keep track of some metrics about our calls
    private CallMetrics callsMetrics;

    // are we being called by another walker (which gets around the traversal-level filtering)?
    private boolean calledByAnotherWalker = true;

    /** Enable deletions in the pileup **/
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /**
     * Initialize the samples, output, and genotype calculation model
     *
     **/
    public void initialize() {

        if ( UAC.POOLED && UAC.genotypeModel == GenotypeCalculationModel.Model.EM_POINT_ESTIMATE ) {
            throw new StingException("This was an attempt to use an EM Point Estimate model with pooled genotype calculations. This model does not work with pooled data.");
        }

        // get all of the unique sample names
        samples = new HashSet<String>();
        List<SAMReadGroupRecord> readGroups = getToolkit().getSAMFileHeader().getReadGroups();
        for ( SAMReadGroupRecord readGroup : readGroups )
            samples.add(readGroup.getSample());

        // print them out for debugging (need separate loop to ensure uniqueness)
        for ( String sample : samples )
            logger.debug("SAMPLE: " + sample);

        gcm = GenotypeCalculationModelFactory.makeGenotypeCalculation(samples, logger, UAC);

        // *** If we were called by another walker, then we don't ***
        // *** want to do any of the other initialization steps.  ***
        if ( VARIANTS_FILE == null && out == null )
            return;

        // if we got here, then we were instantiated by the GATK engine
        calledByAnotherWalker = false;

        // create the output writer stream
        if ( VARIANTS_FILE != null )
            writer = GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), VARIANTS_FILE,
                                                  "UnifiedGenotyper",
                                                  this.getToolkit().getArguments().referenceFile.getName(),
                                                  samples);
        else
            writer = GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), out, "UnifiedGenotyper",
                                                  this.getToolkit().getArguments().referenceFile.getName(),
                                                  samples);
        callsMetrics = new CallMetrics();
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param context contextual information around the locus
     */
    public Pair<List<GenotypeCall>, GenotypeMetaData> map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext context) {
        char ref = Character.toUpperCase(refContext.getBase());
        if ( !BaseUtils.isRegularBase(ref) )
            return null;

        // because other walkers externally call this map method with their own contexts,
        // we need to make sure that the reads are appropriately filtered
        if ( calledByAnotherWalker )
            context = filterAlignmentContext(context);

        // an optimization to speed things up when there is no coverage or when overly covered
        if ( context.getReads().size() == 0 ||
             (UAC.MAX_READS_IN_PILEUP > 0 && context.getReads().size() > UAC.MAX_READS_IN_PILEUP) )
            return null;

        DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, UAC.heterozygosity, DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);
        return gcm.calculateGenotype(tracker, ref, context, priors);
    }

    // filter the given alignment context
    private AlignmentContext filterAlignmentContext(AlignmentContext context) {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        List<SAMRecord> newReads = new ArrayList<SAMRecord>();
        List<Integer> newOffsets = new ArrayList<Integer>();

        for (int i = 0; i < reads.size(); i++) {
            SAMRecord read = reads.get(i);
            if ( read.getMappingQuality() != 0 && read.getReadGroup() != null ) {
                newReads.add(read);
                newOffsets.add(offsets.get(i));
            }
        }

        return new AlignmentContext(context.getLocation(), newReads, newOffsets);                
    }

    /**
     * Determine whether we're at a Hapmap site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a Hapmap site, false if otherwise
     */
    private static boolean isHapmapSite(RefMetaDataTracker tracker) {
        return tracker.getTrackData("hapmap", null) != null;
    }

    /**
     * Determine whether we're at a dbSNP site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a dbSNP site, false if otherwise
     */
    private static boolean isDbSNPSite(RefMetaDataTracker tracker) {
        return rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null)) != null;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Pair<List<GenotypeCall>, GenotypeMetaData> value, Integer sum) {
        if ( value == null || value.first == null )
            return sum;

        callsMetrics.nCalledBases++;

        if ( value.first.size() == 0 )
            return sum;

        // special-case for single-sample using PointEstimate model
        if ( value.second == null ) {
            GenotypeCall call = value.first.get(0);
            if ( UAC.GENOTYPE || call.isVariant(call.getReference()) ) {
                double confidence = (UAC.GENOTYPE ? call.getNegLog10PError() : call.toVariation().getNegLog10PError());
                if ( confidence >= UAC.LOD_THRESHOLD ) {
                    callsMetrics.nConfidentCalls++;
                    writer.addGenotypeCall(call);
                }
            } else {
                callsMetrics.nNonConfidentCalls++;
            }
        }

        // use multi-sample mode if we have multiple samples or the output type allows it
        else if ( writer.supportsMultiSample() || samples.size() > 1 ) {

            // annoying hack to get around Java generics
            ArrayList<Genotype> callList = new ArrayList<Genotype>();
            for ( GenotypeCall call : value.first )
                callList.add(call);

            callsMetrics.nConfidentCalls++;
            writer.addMultiSampleCall(callList, value.second);
        }

        // otherwise, use single sample mode
        else {
            callsMetrics.nConfidentCalls++;
            writer.addGenotypeCall(value.first.get(0));
        }

        return sum + 1;
    }

    /** Close the variant file. */
    public void onTraversalDone(Integer sum) {
        logger.info("Processed " + sum + " loci that are callable for SNPs");
        writer.close();
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
