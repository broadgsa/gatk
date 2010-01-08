package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.gatk.walkers.variantstovcf.VariantsToVCF;

import java.util.*;
import java.io.*;


/**
 * Annotates variant calls with context information.  Users can specify which of the available annotations to use.
 */
//@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=VariationRod.class))
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@Reference(window=@Window(start=-20,stop=20))
public class VariantAnnotator extends RodWalker<Integer, Integer> {
    @Argument(fullName="vcfOutput", shortName="vcf", doc="VCF file to which all variants should be written with annotations", required=true)
    protected File VCF_OUT;
    @Argument(fullName="sampleName", shortName="sample", doc="The sample (NA-ID) corresponding to the variant input (for non-VCF input only)", required=false)
    protected String sampleName = null;
    @Argument(fullName="annotations", shortName="A", doc="Annotation types to apply to variant calls", required=false)
    protected String[] ANNOTATIONS;
    @Argument(fullName="includeExperimentalAnnotations", shortName="exp", doc="Use all possible annotations, including experimental ones", required=false)
    protected Boolean USE_ALL_ANNOTATIONS = false;
    @Argument(fullName="useStandardAnnotations", shortName="standard", doc="Use all standard annotations", required=false)
    protected Boolean USE_STANDARD_ANNOTATIONS = false;
    @Argument(fullName="list", shortName="ls", doc="List the available annotations and exit")
    protected Boolean LIST = false;

    private VCFWriter vcfWriter;
    private VCFHeader vcfHeader;

    private HashMap<String, String> nonVCFsampleName = new HashMap<String, String>();

    private ArrayList<VariantAnnotation> requestedAnnotations;

    // mapping from class name to class
    private static HashMap<String, VariantAnnotation> allAnnotations = null;
    private static HashMap<String, VariantAnnotation> standardAnnotations = null;


    private static void determineAllAnnotations() {
        allAnnotations = new HashMap<String, VariantAnnotation>();
        standardAnnotations = new HashMap<String, VariantAnnotation>();
        List<Class<? extends VariantAnnotation>> annotationClasses = PackageUtils.getClassesImplementingInterface(VariantAnnotation.class);
        for ( Class c : annotationClasses ) {
            try {
                VariantAnnotation annot = (VariantAnnotation) c.newInstance();
                allAnnotations.put(c.getSimpleName().toUpperCase(), annot);
                if ( annot instanceof StandardVariantAnnotation )
                    standardAnnotations.put(c.getSimpleName().toUpperCase(), annot);
            } catch (InstantiationException e) {
                throw new StingException(String.format("Cannot instantiate annotation class '%s': must be concrete class", c.getSimpleName()));
            } catch (IllegalAccessException e) {
                throw new StingException(String.format("Cannot instantiate annotation class '%s': must have no-arg constructor", c.getSimpleName()));
            }
        }
    }

    private void listFiltersAndExit() {
        List<Class<? extends VariantAnnotation>> annotationClasses = PackageUtils.getClassesImplementingInterface(VariantAnnotation.class);
        out.println("\nAvailable annotations:");
        for (int i = 0; i < annotationClasses.size(); i++)
            out.println("\t" + annotationClasses.get(i).getSimpleName());
        out.println();
        System.exit(0);
    }

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        if ( LIST )
            listFiltersAndExit();

        // get the list of all sample names from the various VCF input rods
        TreeSet<String> samples = new TreeSet<String>();
        VCFUtils.getUniquifiedSamplesFromRods(getToolkit(), samples, new HashMap<Pair<String, String>, String>());

        // add the non-VCF sample from the command-line, if applicable
        if ( sampleName != null ) {
            nonVCFsampleName.put(sampleName.toUpperCase(), "variant");
            samples.add(sampleName.toUpperCase());
        }

        // if there are no valid samples, warn the user
        if ( samples.size() == 0 )
            logger.warn("There are no samples input at all; use the --sampleName argument to specify one if desired.");

        determineAllAnnotations();

        if ( USE_STANDARD_ANNOTATIONS ) {
            requestedAnnotations = new ArrayList<VariantAnnotation>(standardAnnotations.values());
        } else if ( USE_ALL_ANNOTATIONS ) {
            requestedAnnotations = new ArrayList<VariantAnnotation>(allAnnotations.values());
        } else {
            requestedAnnotations = new ArrayList<VariantAnnotation>();
            if ( ANNOTATIONS != null ) {
                for ( String requested : ANNOTATIONS ) {

                    VariantAnnotation annot = allAnnotations.get(requested.toUpperCase());
                    if ( annot == null )
                        throw new StingException("Unknown annotation '" + requested + "'.  Issue the '-ls' argument to list available annotations.");

                    requestedAnnotations.add(annot);
                }
            }
        }

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "VariantAnnotator"));
        hInfo.add(new VCFHeaderLine("annotatorReference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.addAll(getVCFAnnotationDescriptions(requestedAnnotations));

        vcfWriter = new VCFWriter(VCF_OUT);
        vcfHeader = new VCFHeader(hInfo, samples);
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

        RODRecordList<ReferenceOrderedDatum> rods = tracker.getTrackData("variant", null);
        // ignore places where we don't have a variant
        if ( rods == null || rods.getRecords().size() == 0 )
            return 0;

        Map<String, String> annotations = new HashMap<String, String>();
        VariationRod variant = (VariationRod)rods.getRecords().get(0);

        // if the reference base is not ambiguous, the variant is a SNP, and it's the appropriate type, we can annotate
        if ( BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 &&
                variant.isBiallelic() &&
                variant.isSNP() ) {
            Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
            if ( stratifiedContexts != null )
                annotations = getAnnotations(tracker, ref, stratifiedContexts, variant, requestedAnnotations);
        }
        writeVCF(tracker, ref, context, variant, annotations);

        return 1;
    }

    // option #1: don't specify annotations to be used: standard annotations are used by default
    public static Set<VCFHeaderLine> getVCFAnnotationDescriptions() {
        if ( standardAnnotations == null )
            determineAllAnnotations();

        TreeSet<VCFHeaderLine> descriptions = new TreeSet<VCFHeaderLine>();
        for ( VariantAnnotation annotation : standardAnnotations.values() )
            descriptions.add(annotation.getDescription());

        return descriptions;
    }

    // option #2: specify that all possible annotations be used
    public static Set<VCFHeaderLine> getAllVCFAnnotationDescriptions() {
        if ( standardAnnotations == null )
            determineAllAnnotations();

        TreeSet<VCFHeaderLine> descriptions = new TreeSet<VCFHeaderLine>();
        for ( VariantAnnotation annotation : allAnnotations.values() )
            descriptions.add(annotation.getDescription());

        return descriptions;
    }

    // option #3: specify the exact annotations to be used
    public static Set<VCFHeaderLine> getVCFAnnotationDescriptions(Collection<VariantAnnotation> annotations) {

        TreeSet<VCFHeaderLine> descriptions = new TreeSet<VCFHeaderLine>();
        for ( VariantAnnotation annotation : annotations )
            descriptions.add(annotation.getDescription());

        return descriptions;
    }

    // option #1: don't specify annotations to be used: standard annotations are used by default
    public static Map<String, String> getAnnotations(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {
        if ( standardAnnotations == null )
            determineAllAnnotations();
        return getAnnotations(tracker, ref, stratifiedContexts, variation, standardAnnotations.values());
    }

    // option #2: specify that all possible annotations be used
    public static Map<String, String> getAllAnnotations(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {
        if ( allAnnotations == null )
            determineAllAnnotations();
        return getAnnotations(tracker, ref, stratifiedContexts, variation, allAnnotations.values());
    }

    // option #3: specify the exact annotations to be used
    public static Map<String, String> getAnnotations(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation, Collection<VariantAnnotation> annotations) {

        HashMap<String, String> results = new HashMap<String, String>();

        for ( VariantAnnotation annotator : annotations) {
            String annot = annotator.annotate(tracker, ref, stratifiedContexts, variation);
            if ( annot != null ) {
                results.put(annotator.getKeyName(), annot);
            }
        }

        return results;
    }


    private void writeVCF(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariationRod variant, Map<String, String> annotations) {
        VCFRecord rec = getVCFRecord(tracker, ref, context, variant);
        if ( rec != null ) {
            rec.addInfoFields(annotations);
            // also, annotate dbsnp id if available and not already there
            if ( rec.getID() == null || rec.getID().equals(".") ) {
                rodDbSNP dbsnp = rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null));
                if ( dbsnp != null )
                    rec.setID(dbsnp.getRS_ID());
            }
            vcfWriter.addRecord(rec);
        }
    }


    private VCFRecord getVCFRecord(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariationRod variant) {
        if ( variant instanceof RodVCF ) {
            return ((RodVCF)variant).mCurrentRecord;
        } else {
            List<VCFGenotypeRecord> gt = new ArrayList<VCFGenotypeRecord>();
            Map<VCFHeader.HEADER_FIELDS, String> map = new HashMap<VCFHeader.HEADER_FIELDS, String>();
            return VariantsToVCF.generateVCFRecord(tracker, ref, context, vcfHeader, gt, map, nonVCFsampleName, out, false, false);
        }
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
