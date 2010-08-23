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

import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;
import java.io.PrintStream;
import java.io.File;
import java.io.FileNotFoundException;


/**
 * A variant caller which unifies the approaches of several disparate callers.  Works for single-sample and
 * multi-sample data.  The user can choose from several different incorporated calculation models.
 */
@Reference(window=@Window(start=-20,stop=20))
@By(DataSource.REFERENCE)
public class UnifiedGenotyper extends LocusWalker<VariantCallContext, UnifiedGenotyper.UGStatistics> implements TreeReducible<UnifiedGenotyper.UGStatistics> {
    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Output(doc="File to which variants should be written",required=false)
    public VCFWriter writer = null;

    @Argument(fullName="variants_out",shortName="varout",doc="Please use --out instead",required=false)
    @Deprecated
    protected String varout;

    @Argument(fullName = "verbose_mode", shortName = "verbose", doc = "File to print all of the annotated and detailed debugging output", required = false)
    protected File verboseFile = null;

    @Argument(fullName = "metrics_file", shortName = "metrics", doc = "File to print any relevant callability metrics output", required = false)
    protected File metricsFile = null;

    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<String>();

    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { "Standard" };

    @Argument(fullName = "NO_HEADER", shortName = "NO_HEADER", doc = "Don't output the usual VCF header tag with the command line. FOR DEBUGGING PURPOSES ONLY. This option is required in order to pass integration tests.", required = false)
    protected Boolean NO_VCF_HEADER_LINE = false;


    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable extended events for indels
    public boolean generateExtendedEvents() { return UAC.genotypeModel == GenotypeCalculationModel.Model.INDELS; }

    private PrintStream verboseWriter;
    private PrintStream metricsWriter;

    /**
     * Inner class for collecting output statistics from the UG
     */
    public static class UGStatistics {
        /** The total number of passes examined -- i.e., the number of map calls */
        long nBasesVisited = 0;

        /** The number of bases that were potentially callable -- i.e., those not at excessive coverage or masked with N */
        long nBasesCallable = 0;

        /** The number of bases called confidently (according to user threshold), either ref or other */
        long nBasesCalledConfidently = 0;

        /** The number of bases for which calls were emitted */
        long nCallsMade = 0;

        /** The total number of extended events encountered */
        long nExtendedEvents = 0;

        double percentCallableOfAll()    { return (100.0 * nBasesCallable) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfAll()      { return (100.0 * nBasesCalledConfidently) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfCallable() { return (100.0 * nBasesCalledConfidently) / (nBasesVisited-nExtendedEvents); }
    }

    /**
     * Initialize the samples, output, and genotype calculation model
     *
     **/
    public void initialize() {
        try {
            if(verboseWriter != null) verboseWriter = new PrintStream(verboseFile);
        }
        catch(FileNotFoundException ex) {
            throw new StingException("Unable to create verbose writer from provided file "+verboseFile,ex);
        }
        try {
            if(metricsWriter != null) metricsWriter = new PrintStream(metricsFile);
        }
        catch(FileNotFoundException ex) {
            throw new StingException("Unable to create metrics writer from provided file "+metricsFile,ex);
        }

        annotationEngine = new VariantAnnotatorEngine(getToolkit(), Arrays.asList(annotationClassesToUse), annotationsToUse);
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, writer, verboseWriter, annotationEngine);

        // initialize the writers
        if ( verboseWriter != null ) {
            StringBuilder header = new StringBuilder("AFINFO\tLOC\tMAF\tF\tNullAFpriors\t");
            for ( byte altAllele : BaseUtils.BASES ) {
                header.append("POfDGivenAFFor" + (char)altAllele + "\t");
                header.append("PosteriorAFFor" + (char)altAllele + "\t");
            }
            verboseWriter.println(header);
        }

        // initialize the header
        writer.writeHeader(new VCFHeader(getHeaderInfo(), UG_engine.samples)) ;
    }

    private Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // annotation (INFO) fields from UnifiedGenotyper
        for ( Map.Entry<String, String> dbSet : UG_engine.dbAnnotations.entrySet() )
            headerInfo.add(new VCFInfoHeaderLine(dbSet.getValue(), 0, VCFHeaderLineType.Flag, (dbSet.getKey().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ? "dbSNP" : dbSet.getValue()) + " Membership"));
        if ( !UAC.NO_SLOD )
            headerInfo.add(new VCFInfoHeaderLine(VCFConstants.STRAND_BIAS_KEY, 1, VCFHeaderLineType.Float, "Strand Bias"));

        // FORMAT and INFO fields
        headerInfo.addAll(VCFUtils.getSupportedHeaderStrings());

        // FILTER fields
        if ( UAC.STANDARD_CONFIDENCE_FOR_EMITTING < UAC.STANDARD_CONFIDENCE_FOR_CALLING ||
             UAC.TRIGGER_CONFIDENCE_FOR_EMITTING < UAC.TRIGGER_CONFIDENCE_FOR_CALLING )
            headerInfo.add(new VCFFilterHeaderLine(UnifiedGenotyperEngine.LOW_QUAL_FILTER_NAME, "Low quality"));

        // all of the arguments from the argument collection
        if ( !NO_VCF_HEADER_LINE ) {
            Set<Object> args = new HashSet<Object>();
            args.add(UAC);
            headerInfo.add(new VCFHeaderLine("UnifiedGenotyper", "\"" + CommandLineUtils.createApproximateCommandLineArgumentString(getToolkit(), args, getClass()) + "\""));
        }
        
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
        lhs.nCallsMade += rhs.nCallsMade;
        return lhs;
    }

    public UGStatistics reduce(VariantCallContext value, UGStatistics sum) {
        // we get a point for reaching reduce
        sum.nBasesVisited++;

        // can't call the locus because of no coverage
        if ( value == null )
            return sum;

        // A call was attempted -- the base was potentially callable
        sum.nBasesCallable++;

        // the base was confidently callable
        sum.nBasesCalledConfidently += value.confidentlyCalled ? 1 : 0;

        // can't make a confident variant call here
        if ( value.vc == null )
            return sum;

        try {
            // we are actually making a call
            sum.nCallsMade++;
            writer.add(value.vc, value.refBase);
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
        logger.info(String.format("Actual calls made                            %d", sum.nCallsMade));

        if ( metricsWriter != null ) {
            metricsWriter.println(String.format("Visited bases                                %d", sum.nBasesVisited));
            metricsWriter.println(String.format("Callable bases                               %d", sum.nBasesCallable));
            metricsWriter.println(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
            metricsWriter.println(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
            metricsWriter.println(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
            metricsWriter.println(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
            metricsWriter.println(String.format("Actual calls made                            %d", sum.nCallsMade));
        }
    }
}
