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
import org.broadinstitute.sting.utils.cmdLine.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.util.*;
import java.io.PrintStream;
import java.io.File;


/**
 * A variant caller which unifies the approaches of several disparate callers.  Works for single-sample,
 * multi-sample, and pooled data.  The user can choose from several different incorporated calculation models.
 */
@Reference(window=@Window(start=-20,stop=20))
@By(DataSource.REFERENCE)
public class UnifiedGenotyper extends LocusWalker<VariantCallContext, UnifiedGenotyper.UGStatistics> implements TreeReducible<UnifiedGenotyper.UGStatistics> {
    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Argument(doc = "File to which variants should be written", required = false)
    public GenotypeWriter writer = null;

    @Argument(fullName = "verbose_mode", shortName = "verbose", doc = "File to print all of the annotated and detailed debugging output", required = false)
    public PrintStream verboseWriter = null;

    @Argument(fullName = "beagle_file", shortName = "beagle", doc = "File to print BEAGLE-specific data for use with imputation", required = false)
    public PrintStream beagleWriter = null;

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;

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
     * Initialize the samples, output, and genotype calculation model
     *
     **/
    public void initialize() {

        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, writer, verboseWriter, beagleWriter);

        // initialize the writers
        if ( verboseWriter != null ) {
            StringBuilder header = new StringBuilder("AFINFO\tLOC\tMAF\tF\tNullAFpriors\t");
            for ( char altAllele : BaseUtils.BASES ) {
                char base = Character.toUpperCase(altAllele);
                header.append("POfDGivenAFFor" + base + "\t");
                header.append("PosteriorAFFor" + base + "\t");
            }
            verboseWriter.println(header);
        }
        if ( beagleWriter != null ) {
            beagleWriter.print("marker alleleA alleleB");
            for ( String sample : UG_engine.samples )
                beagleWriter.print(String.format(" %s %s %s", sample, sample, sample));
            beagleWriter.println();
        }

        // initialize the header
        GenotypeWriterFactory.writeHeader(writer, GenomeAnalysisEngine.instance.getSAMFileHeader(), UG_engine.samples, getHeaderInfo());
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
        if ( UG_engine.annotateDbsnp )
            headerInfo.add(new VCFInfoHeaderLine(VCFRecord.DBSNP_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "dbSNP Membership"));
        if ( UG_engine.annotateHapmap2 )
            headerInfo.add(new VCFInfoHeaderLine(VCFRecord.HAPMAP2_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "HapMap2 Membership"));
        if ( UG_engine.annotateHapmap3 )
            headerInfo.add(new VCFInfoHeaderLine(VCFRecord.HAPMAP3_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "HapMap3 Membership"));
        if ( !UAC.NO_SLOD )
            headerInfo.add(new VCFInfoHeaderLine(VCFRecord.STRAND_BIAS_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Strand Bias"));

        // FORMAT and INFO fields if not in POOLED mode
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
        // also, the list of input bams
        for ( File file : getToolkit().getArguments().samFiles )
            headerInfo.add(new VCFHeaderLine("UG_bam_file_used", file.getName()));

        return headerInfo;
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public VariantCallContext map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        return UG_engine.runGenotyper(tracker, refContext, rawContext);
    }

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
        if ( value.vc == null )
            return sum;

        try {
            writer.addCall(value.vc);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(e.getMessage() + "; this is often caused by using the --assume_single_sample_reads argument with the wrong sample name");
        }

        return sum;
    }

    public void onTraversalDone(UGStatistics sum) {
        logger.info(String.format("Visited bases                                %d", sum.nBasesVisited));
        logger.info(String.format("Callable bases                               %d", sum.nBasesCallable));
        logger.info(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
        logger.info(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
        logger.info(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
        logger.info(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
    }
}
