package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.io.*;


/**
 * Annotates variant calls with context information.  Users can specify which of the available annotations to use.
 */
//@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=VariationRod.class))
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@Reference(window=@Window(start=-50,stop=50))
@By(DataSource.REFERENCE)
public class VariantAnnotator extends LocusWalker<Integer, Integer> {
    @Argument(fullName="vcfOutput", shortName="vcf", doc="VCF file to which all variants should be written with annotations", required=true)
    protected File VCF_OUT;
    @Argument(fullName="sampleName", shortName="sample", doc="The sample (NA-ID) corresponding to the variant input (for non-VCF input only)", required=false)
    protected String sampleName = null;
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected String[] annotationsToUse = {};
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { };
    @Argument(fullName="useAllAnnotations", shortName="all", doc="Use all possible annotations (not for the faint of heart)", required=false)
    protected Boolean USE_ALL_ANNOTATIONS = false;
    @Argument(fullName="list", shortName="ls", doc="List the available annotations and exit")
    protected Boolean LIST = false;

    private VCFWriter vcfWriter;

    private HashMap<String, String> nonVCFsampleName = new HashMap<String, String>();

    private VariantAnnotatorEngine engine;


    private void listAnnotationsAndExit() {
        List<Class<? extends InfoFieldAnnotation>> infoAnnotationClasses = PackageUtils.getClassesImplementingInterface(InfoFieldAnnotation.class);
        out.println("\nAvailable annotations for the VCF INFO field:");
        for (int i = 0; i < infoAnnotationClasses.size(); i++)
            out.println("\t" + infoAnnotationClasses.get(i).getSimpleName());
        out.println();
        List<Class<? extends GenotypeAnnotation>> genotypeAnnotationClasses = PackageUtils.getClassesImplementingInterface(GenotypeAnnotation.class);
        out.println("\nAvailable annotations for the VCF FORMAT field:");
        for (int i = 0; i < genotypeAnnotationClasses.size(); i++)
            out.println("\t" + genotypeAnnotationClasses.get(i).getSimpleName());
        out.println();
        System.exit(0);
    }

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        if ( LIST )
            listAnnotationsAndExit();

        // get the list of all sample names from the various VCF input rods
        TreeSet<String> samples = new TreeSet<String>();
        SampleUtils.getUniquifiedSamplesFromRods(getToolkit(), samples, new HashMap<Pair<String, String>, String>());

        // add the non-VCF sample from the command-line, if applicable
        if ( sampleName != null  ) {
            nonVCFsampleName.put(sampleName.toUpperCase(), "variant");
            samples.add(sampleName.toUpperCase());
        }

        // if there are no valid samples, warn the user
        if ( samples.size() == 0 ) {
            logger.warn("There are no samples input at all; use the --sampleName argument to specify one if desired.");
        }

        if ( USE_ALL_ANNOTATIONS )
            engine = new VariantAnnotatorEngine(getToolkit());
        else
            engine = new VariantAnnotatorEngine(getToolkit(), annotationClassesToUse, annotationsToUse);

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "VariantAnnotator"));
        hInfo.add(new VCFHeaderLine("annotatorReference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.addAll(engine.getVCFAnnotationDescriptions());

        vcfWriter = new VCFWriter(VCF_OUT);
        VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
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

        RODRecordList rods = tracker.getTrackData("variant", null);
        // ignore places where we don't have a variant
        if ( rods == null || rods.size() == 0 )
            return 0;

        ReferenceOrderedDatum variant = rods.get(0);
        VariantContext vc = VariantContextAdaptors.toVariantContext("variant", variant);
        if ( vc == null )
            return 0;

        // if the reference base is not ambiguous, we can annotate
        if ( BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 ) {
            Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
            if ( stratifiedContexts != null ) {
                vc = engine.annotateContext(tracker, ref, stratifiedContexts, vc);
            }
        }

        if ( variant instanceof RodVCF ) {
            RodVCF vcf = (RodVCF)variant;
            vcfWriter.addRecord(VariantContextAdaptors.toVCF(vc, ref.getBase(), Arrays.asList(vcf.getRecord().getGenotypeFormatString().split(VCFRecord.GENOTYPE_FIELD_SEPERATOR)), vcf.getFilterString() != null));
        } else {
            vcfWriter.addRecord(VariantContextAdaptors.toVCF(vc, ref.getBase()));
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
        out.printf("Processed %d loci.\n", result);

        vcfWriter.close();
    }

    public static Genotype getFirstVariant(char ref, List<Genotype> genotypes) {
        for ( Genotype g : genotypes ) {
            if ( g.isVariant(ref) )
                return g;
        }
        return null;
    }
}
