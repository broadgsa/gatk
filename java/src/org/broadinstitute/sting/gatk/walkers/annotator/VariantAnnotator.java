package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
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
@Reference(window=@Window(start=-20,stop=20))
@By(DataSource.REFERENCE)
public class VariantAnnotator extends LocusWalker<Integer, Integer> {
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

    private HashMap<String, String> nonVCFsampleName = new HashMap<String, String>();

    private ArrayList<VariantAnnotation> requestedAnnotations;

    // should we annotate dbsnp?
    private boolean annotateDbsnp = false;
    // how about hapmap2?
    private boolean annotateHapmap2 = false;
    // how about hapmap3?
    private boolean annotateHapmap3 = false;

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

        // check to see whether a dbsnp rod was included
        List<ReferenceOrderedDataSource> dataSources = getToolkit().getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            ReferenceOrderedData rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(rodDbSNP.class) ) {
                annotateDbsnp = true;
            }
            if ( rod.getName().equals("hapmap2") ) {
                annotateHapmap2 = true;
            }
            if ( rod.getName().equals("hapmap3") ) {
                annotateHapmap3 = true;
            }
        }

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "VariantAnnotator"));
        hInfo.add(new VCFHeaderLine("annotatorReference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.addAll(getVCFAnnotationDescriptions(requestedAnnotations));
        if ( annotateDbsnp )
            hInfo.add(new VCFInfoHeaderLine(VCFRecord.DBSNP_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "dbSNP membership"));
        if ( annotateHapmap2 )
            hInfo.add(new VCFInfoHeaderLine(VCFRecord.HAPMAP2_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Hapmap 2 membership"));
        if ( annotateHapmap3 )
            hInfo.add(new VCFInfoHeaderLine(VCFRecord.HAPMAP3_KEY,1,VCFInfoHeaderLine.INFO_TYPE.Integer, "Hapmap 3 membership"));

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

        Map<String, String> annotations = new HashMap<String, String>();
        ReferenceOrderedDatum variant = rods.get(0);
        VariantContext vc = VariantContextAdaptors.toVariantContext("variant", variant);
        if ( vc == null )
            return 0;

        // if the reference base is not ambiguous, we can annotate
        if ( BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 ) {
            Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
            if ( stratifiedContexts != null )
                annotations = getAnnotations(tracker, ref, stratifiedContexts, vc, requestedAnnotations, annotateDbsnp, annotateHapmap2, annotateHapmap3);
        }

        VCFRecord record;
        if ( variant instanceof RodVCF )
            record = ((RodVCF)variant).mCurrentRecord;
        else
            record = VariantContextAdaptors.toVCF(vc, ref.getBase());

        record.addInfoFields(annotations);
        writeVCF(tracker, record);

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
    public static Map<String, String> getAnnotations(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc, boolean annotateDbsnp, boolean annotateHapmap2, boolean annotateHapmap3) {
        if ( standardAnnotations == null )
            determineAllAnnotations();
        return getAnnotations(tracker, ref, stratifiedContexts, vc, standardAnnotations.values(), annotateDbsnp, annotateHapmap2, annotateHapmap3);
    }

    // option #2: specify that all possible annotations be used
    public static Map<String, String> getAllAnnotations(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc, boolean annotateDbsnp, boolean annotateHapmap2, boolean annotateHapmap3) {
        if ( allAnnotations == null )
            determineAllAnnotations();
        return getAnnotations(tracker, ref, stratifiedContexts, vc, allAnnotations.values(), annotateDbsnp, annotateHapmap2, annotateHapmap3);
    }

    // option #3: specify the exact annotations to be used
    public static Map<String, String> getAnnotations(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc, Collection<VariantAnnotation> annotations, boolean annotateDbsnp, boolean annotateHapmap2, boolean annotateHapmap3) {

        HashMap<String, String> results = new HashMap<String, String>();

        // annotate dbsnp occurrence
        if ( annotateDbsnp ) {
            rodDbSNP dbsnp = rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null));
            results.put(VCFRecord.DBSNP_KEY, dbsnp == null ? "0" : "1");
        }

        if ( annotateHapmap2 ) {
            RODRecordList hapmap2 = tracker.getTrackData("hapmap2",null);
            results.put(VCFRecord.HAPMAP2_KEY, hapmap2 == null? "0" : "1");
        }

        if ( annotateHapmap3 ) {
            RODRecordList hapmap3 = tracker.getTrackData("hapmap3",null);
            results.put( VCFRecord.HAPMAP3_KEY, hapmap3 == null ? "0" : "1");
        }

        for ( VariantAnnotation annotator : annotations) {
            String annot = annotator.annotate(tracker, ref, stratifiedContexts, vc);
            if ( annot != null ) {
                results.put(annotator.getKeyName(), annot);
            }
        }

        return results;
    }

    private void writeVCF(RefMetaDataTracker tracker, VCFRecord record) {
        // annotate dbsnp id if available and not already there
        if ( annotateDbsnp && (record.getID() == null || record.getID().equals(VCFRecord.EMPTY_ID_FIELD)) ) {
            rodDbSNP dbsnp = rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null));
            if ( dbsnp != null )
                record.setID(dbsnp.getRS_ID());
        }
        vcfWriter.addRecord(record);
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
