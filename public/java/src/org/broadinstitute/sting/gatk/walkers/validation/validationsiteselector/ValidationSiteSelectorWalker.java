/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package org.broadinstitute.sting.gatk.walkers.validation.validationsiteselector;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.codecs.vcf.writer.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.util.*;


/**
 * Randomly selects VCF records according to specified options.
 *
 * <p>
 * ValidationSiteSelectorWalker is intended for use in experiments where we sample data randomly from a set of variants, for example
 * in order to choose sites for a follow-up validation study.
 *
 * Sites are selected randomly but within certain restrictions. There are two main sources of restrictions
 * a) Sample restrictions. A user can specify a set of samples, and we will only consider sites which are polymorphic within such given sample subset.
 * These sample restrictions can be given as a set of individual samples, a text file (each line containing a sample name), or a regular expression.
 * A user can additionally specify whether samples will be considered based on their genotypes (a non-reference genotype means that such sample is polymorphic in that variant,
 * and hence that variant will be considered for inclusion in set), or based on their PLs.
 * b) A user can additionally specify a sampling method based on allele frequency. Two sampling methods are currently supported.
 * 1. Uniform sampling will just sample uniformly from variants polymorphic in selected samples.
 * 2. Sampling based on Allele Frequency spectrum will ensure that output sites have the same AF distribution as the input set.
 *
 * User can additionally restrict output to a particular type of variant (SNP, Indel, etc.)
 *
 * <h2>Input</h2>
 * <p>
 * One or more variant sets to choose from.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A sites-only VCF with the desired number of randomly selected sites.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ValidationSiteSelectorWalker \
 *   --variant input1.vcf \
 *   --variant input2.vcf \
 *   -sn NA12878 \
 *   -o output.vcf \
 *   --numValidationSites 200   \
 *   -sampleMode  POLY_BASED_ON_GT \
 *   -freqMode KEEP_AF_SPECTRUM
 *
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ValidationSiteSelectorWalker \
 *   --variant:foo input1.vcf \
 *   --variant:bar input2.vcf \
 *   --numValidationSites 200 \
 *   -sf samples.txt \
 *   -o output.vcf \
 *   -sampleMode  POLY_BASED_ON_GT \
  *   -freqMode UNIFORM
 *   -selectType INDEL
 * </pre>
 *
 */
public class ValidationSiteSelectorWalker extends RodWalker<Integer, Integer> {

    public enum AF_COMPUTATION_MODE {
        KEEP_AF_SPECTRUM,
        UNIFORM
    }

    public enum SAMPLE_SELECTION_MODE {
        NONE,
        POLY_BASED_ON_GT,
        POLY_BASED_ON_GL
    }

    /**
     * The input VCF file
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file, can be specified multiple times", required=true)
    public List<RodBinding<VariantContext>> variants;

    /**
     * The output VCF file
     */
    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    /**
     * Sample name(s) to subset the input VCF to, prior to selecting variants. -sn A -sn B subsets to samples A and B.
     */
    @Argument(fullName="sample_name", shortName="sn", doc="Include genotypes from this sample. Can be specified multiple times", required=false)
    public Set<String> sampleNames = new HashSet<String>(0);

    /**
     * Sample regexps to subset the input VCF to, prior to selecting variants. -sn NA12* subsets to all samples with prefix NA12
     */
    @Argument(fullName="sample_expressions", shortName="se", doc="Regular expression to select many samples from the ROD tracks provided. Can be specified multiple times", required=false)
    public Set<String> sampleExpressions ;

    /**
     * File containing a list of sample names to subset the input vcf to. Equivalent to specifying the contents of the file separately with -sn
     */
    @Input(fullName="sample_file", shortName="sf", doc="File containing a list of samples (one per line) to include. Can be specified multiple times", required=false)
    public Set<File> sampleFiles;

    /**
     * A mode for selecting sites based on sample-level data. See the wiki documentation for more information.
     */
    @Argument(fullName="sampleMode", shortName="sampleMode", doc="Sample selection mode", required=false)
    private SAMPLE_SELECTION_MODE sampleMode = SAMPLE_SELECTION_MODE.NONE;

    /**
     * An P[nonref] threshold for SAMPLE_SELECTION_MODE=POLY_BASED_ON_GL. See the wiki documentation for more information.
     */
    @Argument(shortName="samplePNonref",fullName="samplePNonref", doc="GL-based selection mode only: the probability" +
            " that a site is non-reference in the samples for which to include the site",required=false)
    private double samplePNonref = 0.99;

    /**
     * The number of sites in your validation set
     */
    @Argument(fullName="numValidationSites", shortName="numSites", doc="Number of output validation sites", required=true)
    private int numValidationSites;

    /**
     * Do not exclude filtered sites (e.g. not PASS or .) from consideration for validation
     */
    @Argument(fullName="includeFilteredSites", shortName="ifs", doc="If true, will include filtered sites in set to choose variants from", required=false)
    private boolean INCLUDE_FILTERED_SITES = false;

    /**
     * Argument for the frequency selection mode. (AC/AF/AN) are taken from VCF info field, not recalculated. Typically specified for sites-only VCFs that still have AC/AF/AN information.
     */
    @Argument(fullName="ignoreGenotypes", shortName="ignoreGenotypes", doc="If true, will ignore genotypes in VCF, will take AC,AF from annotations and will make no sample selection", required=false)
    private boolean IGNORE_GENOTYPES = false;

    /**
     * Argument for the frequency selection mode. Allows reference (non-polymorphic) sites to be included in the validation set.
     */
    @Argument(fullName="ignorePolymorphicStatus", shortName="ignorePolymorphicStatus", doc="If true, will ignore polymorphic status in VCF, and will take VCF record directly without pre-selection", required=false)
    private boolean IGNORE_POLYMORPHIC = false;

    @Hidden
    @Argument(fullName="numFrequencyBins", shortName="numBins", doc="Number of frequency bins if we're to match AF distribution", required=false)
    private int numFrequencyBins = 20;

    /**
      * This argument selects allele frequency selection mode. See the wiki for more information.
      */
    @Argument(fullName="frequencySelectionMode", shortName="freqMode", doc="Allele Frequency selection mode", required=false)
    private AF_COMPUTATION_MODE freqMode = AF_COMPUTATION_MODE.KEEP_AF_SPECTRUM;

    /**
      * This argument selects particular kinds of variants (i.e. SNP, INDEL) out of a list. If left unspecified, all types are considered.
      */
     @Argument(fullName="selectTypeToInclude", shortName="selectType", doc="Select only a certain type of variants from the input file. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION. Can be specified multiple times", required=false)
     private List<VariantContext.Type> TYPES_TO_INCLUDE = new ArrayList<VariantContext.Type>();


    private TreeSet<String> samples = new TreeSet<String>();
    SampleSelector sampleSelector = null;
    FrequencyModeSelector frequencyModeSelector = null;
    private ArrayList<VariantContext.Type> selectedTypes = new ArrayList<VariantContext.Type>();

    public void initialize() {
         // Get list of samples to include in the output
         Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit());
         TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));

         Collection<String> samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFiles);
         Collection<String> samplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, sampleExpressions);

         // first, add any requested samples
         samples.addAll(samplesFromFile);
         samples.addAll(samplesFromExpressions);
         samples.addAll(sampleNames);

         // if none were requested, we want all of them
         if ( samples.isEmpty() ) {
             samples.addAll(vcfSamples);

         }

         sampleSelector = getSampleSelectorObject(sampleMode, samples);

        // initialize frequency mode selector
        frequencyModeSelector = getFrequencyModeSelectorObject(freqMode, getToolkit().getGenomeLocParser());

        // if user specified types to include, add these, otherwise, add all possible variant context types to list of vc types to include
        if (TYPES_TO_INCLUDE.isEmpty()) {

            for (VariantContext.Type t : VariantContext.Type.values())
                selectedTypes.add(t);

        }
        else {
            for (VariantContext.Type t : TYPES_TO_INCLUDE)
                selectedTypes.add(t);

        }

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.add(new VCFHeaderLine("source", "ValidationSiteSelector"));
        vcfWriter.writeHeader(new VCFHeader(headerLines));

    }


    @Override
     public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
         if ( tracker == null )
             return 0;

        Collection<VariantContext> vcs = tracker.getValues(variants, context.getLocation());

         if ( vcs == null || vcs.size() == 0) {
             return 0;
         }


        for (VariantContext vc : vcs) {
            if (!selectedTypes.contains(vc.getType()))
                continue;

            // skip if site isn't polymorphic and if user didn't request to ignore polymorphic status
            if (!vc.isPolymorphicInSamples() && !IGNORE_POLYMORPHIC)
                continue;

            if (!INCLUDE_FILTERED_SITES && vc.filtersWereApplied() && vc.isFiltered())
                continue;


            // does this site pass the criteria for the samples we are interested in?
            boolean passesSampleSelectionCriteria;
            if (samples.isEmpty())
                passesSampleSelectionCriteria = true;
            else
                passesSampleSelectionCriteria = sampleSelector.selectSiteInSamples(vc);

            frequencyModeSelector.logCurrentSiteData(vc,passesSampleSelectionCriteria,IGNORE_GENOTYPES,IGNORE_POLYMORPHIC);
        }
        return 1;
    }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    public void onTraversalDone(Integer result) {
        logger.info("Outputting validation sites...");
        ArrayList<VariantContext> selectedSites = frequencyModeSelector.selectValidationSites(numValidationSites);

        for (VariantContext vc : selectedSites) {
            vcfWriter.add(vc);
        }
        logger.info(result + " records processed.");

    }

    private SampleSelector getSampleSelectorObject(SAMPLE_SELECTION_MODE sampleMode, TreeSet<String> samples) {
        SampleSelector sm;
         switch ( sampleMode ) {
             case POLY_BASED_ON_GL:
                 sm = new GLBasedSampleSelector(samples, Math.log10(1.0-samplePNonref));
                 break;
             case POLY_BASED_ON_GT:
                 sm = new GTBasedSampleSelector(samples);
                 break;
             case NONE:
                 sm = new NullSampleSelector(samples);
                 break;
             default:
                 throw new IllegalArgumentException("Unsupported Sample Selection Mode: " + sampleMode);
         }

         return sm;
    }

    private FrequencyModeSelector getFrequencyModeSelectorObject (AF_COMPUTATION_MODE freqMode, GenomeLocParser parser) {
        FrequencyModeSelector fm;

        switch (freqMode) {
            case KEEP_AF_SPECTRUM:
                fm = new KeepAFSpectrumFrequencySelector(numFrequencyBins, parser);
                break;
            case UNIFORM:
                fm = new UniformSamplingFrequencySelector(parser);
                break;
            default: throw new IllegalArgumentException("Unexpected Frequency Selection Mode: "+ freqMode);

        }
        return fm;
    }
}
