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

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.*;

/**
 * A variant caller which unifies the approaches of several disparate callers -- Works for single-sample and multi-sample data.
 *
 * <p>
 * The GATK Unified Genotyper is a multiple-sample, technology-aware SNP and indel caller. It uses a Bayesian genotype
 * likelihood model to estimate simultaneously the most likely genotypes and allele frequency in a population of N samples,
 * emitting an accurate posterior probability of there being a segregating variant allele at each locus as well as for the
 * genotype of each sample. The system can either emit just the variant sites or complete genotypes (which includes
 * homozygous reference calls) satisfying some phred-scaled confidence value. The genotyper can make accurate calls on
 * both single sample data and multi-sample data.
 *
 * <h2>Input</h2>
 * <p>
 * The read data from which to make variant calls.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A raw, unfiltered, highly specific callset in VCF format.
 * </p>
 *
 * <h2>Example generic command for multi-sample SNP calling</h2>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -T UnifiedGenotyper \
 *   -I sample1.bam [-I sample2.bam ...] \
 *   --dbsnp dbSNP.vcf \
 *   -o snps.raw.vcf \
 *   -stand_call_conf [50.0] \
 *   -stand_emit_conf 10.0 \
 *   -dcov [50] \
 *   [-L targets.interval_list]
 * </pre>
 *
 * <p>
 * The above command will call all of the samples in your provided BAM files [-I arguments] together and produce a VCF file
 * with sites and genotypes for all samples. The easiest way to get the dbSNP file is from the GATK resource bundle. Several
 * arguments have parameters that should be chosen based on the average coverage per sample in your data. See the detailed
 * argument descriptions below.
 * </p>
 *
 * <h2>Example command for generating calls at all sites</h2>
 * <pre>
 * java -jar /path/to/GenomeAnalysisTK.jar \
 *   -l INFO \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -T UnifiedGenotyper \
 *   -I /DCC/ftp/pilot_data/data/NA12878/alignment/NA12878.SLX.maq.SRP000031.2009_08.bam \
 *   -o my.vcf \
 *   --output_mode EMIT_ALL_SITES
 * </pre>
 *
 * <h2>Caveats</h2>
 * <ul>
 * <li>The system is under active and continuous development. All outputs, the underlying likelihood model, arguments, and
 * file formats are likely to change.</li>
 * <li>The system can be very aggressive in calling variants. In the 1000 genomes project for pilot 2 (deep coverage of ~35x)
 * we expect the raw Qscore > 50 variants to contain at least ~10% FP calls. We use extensive post-calling filters to eliminate
 * most of these FPs. Variant Quality Score Recalibration is a tool to perform this filtering.</li>
 * <li>We only handle diploid genotypes</li>
 * </ul>
 *
 */

@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@ReadFilters( {BadMateFilter.class, MappingQualityUnavailableFilter.class} )
@Reference(window=@Window(start=-200,stop=200))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.BY_SAMPLE, toCoverage=250)
public class UnifiedGenotyper extends LocusWalker<VariantCallContext, UnifiedGenotyper.UGStatistics> implements TreeReducible<UnifiedGenotyper.UGStatistics>, AnnotatorCompatibleWalker {

    @ArgumentCollection
    private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    /**
     * rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }
    public RodBinding<VariantContext> getVariantRodBinding() { return null; }
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getCompRodBindings() { return Collections.emptyList(); }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }

    /**
     * A raw, unfiltered, highly specific callset in VCF format.
     */
    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter writer = null;

    @Argument(fullName = "debug_file", shortName = "debug_file", doc = "File to print all of the annotated and detailed debugging output", required = false)
    protected PrintStream verboseWriter = null;

    @Argument(fullName = "metrics_file", shortName = "metrics", doc = "File to print any relevant callability metrics output", required = false)
    protected PrintStream metricsWriter = null;

    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     */
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<String>();

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<String>();

    /**
     * Which groups of annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available groups.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { "Standard" };

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable extended events for indels
    public boolean generateExtendedEvents() {
        return (UAC.GLmodel != GenotypeLikelihoodsCalculationModel.Model.SNP && UAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES);
    }

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
        double percentCalledOfCallable() { return (100.0 * nBasesCalledConfidently) / (nBasesCallable); }
    }

    /**
     * Initialize the samples, output, and genotype calculation model
     *
     **/
    public void initialize() {
        // get all of the unique sample names
        // if we're supposed to assume a single sample, do so
        Set<String> samples = new TreeSet<String>();
        if ( UAC.ASSUME_SINGLE_SAMPLE != null )
            samples.add(UAC.ASSUME_SINGLE_SAMPLE);
        else
            samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());

        // initialize the verbose writer
        if ( verboseWriter != null )
            verboseWriter.println("AFINFO\tLOC\tREF\tALT\tMAF\tF\tAFprior\tAFposterior\tNormalizedPosterior");

        annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, verboseWriter, annotationEngine, samples);

        // initialize the header
        Set<VCFHeaderLine> headerInfo = getHeaderInfo();

        // invoke initialize() method on each of the annotation classes, allowing them to add their own header lines
        // and perform any necessary initialization/validation steps
        annotationEngine.invokeAnnotationInitializationMethods(headerInfo);

        writer.writeHeader(new VCFHeader(headerInfo, samples));
    }

    private Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // annotation (INFO) fields from UnifiedGenotyper
        if ( UAC.COMPUTE_SLOD )
            headerInfo.add(new VCFInfoHeaderLine(VCFConstants.STRAND_BIAS_KEY, 1, VCFHeaderLineType.Float, "Strand Bias"));
        headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DOWNSAMPLED_KEY, 0, VCFHeaderLineType.Flag, "Were any of the samples downsampled?"));

        // also, check to see whether comp rods were included
        if ( dbsnp.dbsnp.isBound() )
            headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DBSNP_KEY, 0, VCFHeaderLineType.Flag, "dbSNP Membership"));

        // FORMAT and INFO fields
        headerInfo.addAll(getSupportedHeaderStrings());

        // FILTER fields
        if ( UAC.STANDARD_CONFIDENCE_FOR_EMITTING < UAC.STANDARD_CONFIDENCE_FOR_CALLING )
            headerInfo.add(new VCFFilterHeaderLine(UnifiedGenotyperEngine.LOW_QUAL_FILTER_NAME, "Low quality"));

        return headerInfo;
    }

    /**
     * return a set of supported format lines; what we currently support for output in the genotype fields of a VCF
     * @return a set of VCF format lines
     */
    private static Set<VCFFormatHeaderLine> getSupportedHeaderStrings() {
        Set<VCFFormatHeaderLine> result = new HashSet<VCFFormatHeaderLine>();
        result.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
        result.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Float, "Genotype Quality"));
        result.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Read Depth (only filtered reads used for calling)"));
        result.add(new VCFFormatHeaderLine(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"));

        return result;
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
        return UG_engine.calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext);
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

        // can't make a call here
        if ( !value.shouldEmit )
            return sum;

        try {
            // we are actually making a call
            sum.nCallsMade++;
            writer.add(value);
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
