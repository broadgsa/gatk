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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotationType;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;


/**
 * Annotates variant calls with context information.
 *
 * <p>
 * VariantAnnotator is a GATK tool for annotating variant calls based on their context.
 * The tool is modular; new annotations can be written easily without modifying VariantAnnotator itself.
 *
 * <h2>Input</h2>
 * <p>
 * A variant set to annotate and optionally one or more BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * An annotated VCF.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T VariantAnnotator \
 *   -I input.bam \
 *   -o output.vcf \
 *   -A DepthOfCoverage
 *   --variant input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@Requires(value={})
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@Reference(window=@Window(start=-50,stop=50))
@By(DataSource.REFERENCE)
public class VariantAnnotator extends RodWalker<Integer, Integer> implements AnnotatorCompatibleWalker {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();
    public RodBinding<VariantContext> getVariantRodBinding() { return variantCollection.variants; }

    /**
     * The INFO field will be annotated with information on the most biologically-significant effect
     * listed in the SnpEff output file for each variant.
     */
    @Input(fullName="snpEffFile", shortName = "snpEffFile", doc="A SnpEff output file from which to add annotations", required=false)
    public RodBinding<VariantContext> snpEffFile;
    public RodBinding<VariantContext> getSnpEffRodBinding() { return snpEffFile; }

    /**
      * rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.
      */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    /**
      * If a record in the 'variant' track overlaps with a record from the provided comp track, the INFO field will be annotated
      *  as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field).  Records that are filtered in the comp track will be ignored.
      *  Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
      */
    @Input(fullName="comp", shortName = "comp", doc="comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }

    /**
      * An external resource VCF file or files from which to annotate.
      *
      * One can add annotations from one of the resource VCFs to the output.
      * For example, if you want to annotate your 'variant' VCF with the AC field value from the rod bound to 'resource',
      * you can specify '-E resource.AC' and records in the output VCF will be annotated with 'resource.AC=N' when a record exists in that rod at the given position.
      * If multiple records in the rod overlap the given position, one is chosen arbitrarily.
      */
    @Input(fullName="resource", shortName = "resource", doc="external resource VCF file", required=false)
    public List<RodBinding<VariantContext>> resources = Collections.emptyList();
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return resources; }

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    /**
     * See the -list argument to view available annotations.
     */
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<String>();

    /**
     * Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<String>();

    /**
     * See the -list argument to view available groups.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected List<String> annotationGroupsToUse = new ArrayList<String>();

    /**
     * This option enables you to add annotations from one VCF to another.
     *
     * For example, if you want to annotate your 'variant' VCF with the AC field value from the rod bound to 'resource',
     * you can specify '-E resource.AC' and records in the output VCF will be annotated with 'resource.AC=N' when a record exists in that rod at the given position.
     * If multiple records in the rod overlap the given position, one is chosen arbitrarily.
     */
    @Argument(fullName="expression", shortName="E", doc="One or more specific expressions to apply to variant calls; see documentation for more details", required=false)
    protected List<String> expressionsToUse = new ArrayList<String>();

    /**
     * Note that the -XL argument can be used along with this one to exclude annotations.
     */
    @Argument(fullName="useAllAnnotations", shortName="all", doc="Use all possible annotations (not for the faint of heart)", required=false)
    protected Boolean USE_ALL_ANNOTATIONS = false;

    @Argument(fullName="list", shortName="ls", doc="List the available annotations and exit")
    protected Boolean LIST = false;

    @Hidden
    @Argument(fullName = "assume_single_sample_reads", shortName = "single_sample", doc = "The single sample that we should assume is represented in the input bam (and therefore associate with all reads regardless of whether they have read groups)", required = false)
    protected String ASSUME_SINGLE_SAMPLE = null;

    @Hidden
    @Argument(fullName="vcfContainsOnlyIndels", shortName="dels",doc="Use if you are annotating an indel vcf, currently VERY experimental", required = false)
    protected boolean indelsOnly = false;

    @Argument(fullName="family_string",shortName="family",required=false,doc="A family string of the form mom+dad=child for use with the mendelian violation ratio annotation")
    public String familyStr = null;

    @Argument(fullName="MendelViolationGenotypeQualityThreshold",shortName="mvq",required=false,doc="The genotype quality treshold in order to annotate mendelian violation ratio")
    public double minGenotypeQualityP = 0.0;

    private VariantAnnotatorEngine engine;

    private Collection<VariantContext> indelBufferContext;


    private void listAnnotationsAndExit() {
        List<Class<? extends InfoFieldAnnotation>> infoAnnotationClasses = new PluginManager<InfoFieldAnnotation>(InfoFieldAnnotation.class).getPlugins();
        System.out.println("\nAvailable annotations for the VCF INFO field:");
        for (int i = 0; i < infoAnnotationClasses.size(); i++)
            System.out.println("\t" + infoAnnotationClasses.get(i).getSimpleName());
        System.out.println();
        List<Class<? extends GenotypeAnnotation>> genotypeAnnotationClasses = new PluginManager<GenotypeAnnotation>(GenotypeAnnotation.class).getPlugins();
        System.out.println("\nAvailable annotations for the VCF FORMAT field:");
        for (int i = 0; i < genotypeAnnotationClasses.size(); i++)
            System.out.println("\t" + genotypeAnnotationClasses.get(i).getSimpleName());
        System.out.println();
        System.out.println("\nAvailable classes/groups of annotations:");
        for ( Class c : new PluginManager<AnnotationType>(AnnotationType.class).getInterfaces() )
            System.out.println("\t" + c.getSimpleName());
        System.out.println();
        System.exit(0);
    }

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        if ( LIST )
            listAnnotationsAndExit();

        // get the list of all sample names from the variant VCF input rod, if applicable
        List<String> rodName = Arrays.asList(variantCollection.variants.getName());
        Set<String> samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), rodName);

        // if there are no valid samples, warn the user
        if ( samples.size() == 0 ) {
            logger.warn("There are no samples input at all; use the --sampleName argument to specify one if desired.");
        }

        if ( USE_ALL_ANNOTATIONS )
            engine = new VariantAnnotatorEngine(annotationsToExclude, this, getToolkit());
        else
            engine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse, annotationsToExclude, this, getToolkit());
        engine.initializeExpressions(expressionsToUse);

        // setup the header fields
        // note that if any of the definitions conflict with our new ones, then we want to overwrite the old ones
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(engine.getVCFAnnotationDescriptions());
        for ( VCFHeaderLine line : VCFUtils.getHeaderFields(getToolkit(), Arrays.asList(variantCollection.variants.getName())) ) {
            if ( isUniqueHeaderLine(line, hInfo) )
                hInfo.add(line);
        }
        for ( String expression : expressionsToUse )
            hInfo.add(new VCFInfoHeaderLine(expression, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Value transferred from another external VCF resource"));

        engine.invokeAnnotationInitializationMethods(hInfo);

        VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);

        if ( indelsOnly ) {
            indelBufferContext = null;
        }
    }

    public static boolean isUniqueHeaderLine(VCFHeaderLine line, Set<VCFHeaderLine> currentSet) {
        if ( !(line instanceof VCFCompoundHeaderLine) )
            return true;

        for ( VCFHeaderLine hLine : currentSet ) {
            if ( hLine instanceof VCFCompoundHeaderLine && ((VCFCompoundHeaderLine)line).sameLineTypeAndName((VCFCompoundHeaderLine)hLine) )
                return false;
        }

        return true;
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }


    /**
     * We want reads that span deletions
     *
     * @return true
     */
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /**
     * We want to see extended events if annotating indels
     *
     * @return true
     */
    public boolean generateExtendedEvents() { return indelsOnly; }

    /**
     * For each site of interest, annotate based on the requested annotation types
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
        if ( VCs.size() == 0 )
            return 0;

        Collection<VariantContext> annotatedVCs = VCs;

        // if the reference base is not ambiguous, we can annotate
        Map<String, AlignmentContext> stratifiedContexts;
        if ( BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 ) {
            if ( ! context.hasExtendedEventPileup() ) {
                stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(context.getBasePileup(), ASSUME_SINGLE_SAMPLE);
            } else {
                stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(context.getExtendedEventPileup(), ASSUME_SINGLE_SAMPLE);
            }
            if ( stratifiedContexts != null ) {
                annotatedVCs = new ArrayList<VariantContext>(VCs.size());
                for ( VariantContext vc : VCs )
                    annotatedVCs.add(engine.annotateContext(tracker, ref, stratifiedContexts, vc));
            }
        }

        if ( ! indelsOnly ) {
            for ( VariantContext annotatedVC : annotatedVCs )
                vcfWriter.add(annotatedVC);
        } else {
            // check to see if the buffered context is different (in location) this context
            if ( indelBufferContext != null && ! VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(),indelBufferContext.iterator().next()).equals(VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(),annotatedVCs.iterator().next())) ) {
                for ( VariantContext annotatedVC : indelBufferContext )
                    vcfWriter.add(annotatedVC);
                indelBufferContext = annotatedVCs;
            } else {
                indelBufferContext = annotatedVCs;
            }
        }

        return 1;
    }

    /**
     * Increment the number of loci processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of loci processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     *
     * @param result  the number of loci seen.
     */
    public void onTraversalDone(Integer result) {
        logger.info("Processed " + result + " loci.\n");
    }
}
