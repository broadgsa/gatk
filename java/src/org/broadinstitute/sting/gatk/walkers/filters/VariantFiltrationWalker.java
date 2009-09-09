package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.geli.GeliTextWriter;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.gatk.walkers.variantstovcf.VariantsToVCF;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.io.*;


/**
 * VariantFiltrationWalker applies specified conditionally independent features and filters to pre-called variants.
 * The former modifiesthe likelihoods of each genotype, while the latter makes a decision to include or exclude a
 * variant outright.  At the moment, the variants are expected to be in gelitext format.
 */
@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type= RodGeliText.class))
public class VariantFiltrationWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="vcfOutput", shortName="vcf", doc="VCF file to which all variants should be written with annotations", required=true) public File VCF_OUT;
    @Argument(fullName="sampleName", shortName="sample", doc="Temporary hack to get VCF to work: the sample (NA-ID) corresponding to the variants", required=true) public String sampleName;
    @Argument(fullName="includedOutput", shortName="included", doc="File to which all variants passing filters should be written", required=true) public String INCLUDED_OUT;
    @Argument(fullName="annotatedOutput", shortName="annotated", doc="File to which all variants should be written with annotations - for debugging/parameterizing", required=false) public String ANNOTATED_OUT;
    @Argument(fullName="features", shortName="F", doc="Feature test (optionally with arguments) to apply to genotype posteriors.  Syntax: 'testname[:arguments]'", required=false) public String[] FEATURES;
    @Argument(fullName="exclusion_criterion", shortName="X", doc="Exclusion test (optionally with arguments) to apply to variant call.  Syntax: 'testname[:key1=arg1,key2=arg2,...]'", required=false) public String[] EXCLUSIONS;
    @Argument(fullName="inclusion_threshold", shortName="IT", doc="The product of the probability to include variants based on these filters must be greater than the value specified here in order to be included", required=false) public Double INCLUSION_THRESHOLD = 0.9;
    @Argument(fullName="verbose", shortName="V", doc="Show how the variant likelihoods are changing with the application of each feature") public Boolean VERBOSE = false;
    @Argument(fullName="list", shortName="ls", doc="List the available features and exclusion criteria and exit") public Boolean LIST = false;

    private List<Class<? extends IndependentVariantFeature>> featureClasses;
    private List<Class<? extends VariantExclusionCriterion>> exclusionClasses;

    private VCFWriter vcfWriter;
    private VCFHeader vcfHeader;
    private PrintWriter annotatedWriter;
    private PrintWriter includedWriter;
    private HashMap<String, String> sampleNames = new HashMap<String, String>();

    private ArrayList<IndependentVariantFeature> requestedFeatures;
    private ArrayList<VariantExclusionCriterion> requestedExclusions;

    // the structures necessary to initialize and maintain a windowed context
    private VariantContextWindow variantContextWindow;
    private static final int windowSize = 10;  // 10 variants on either end of the current one
    private ArrayList<VariantContext> windowInitializer = new ArrayList<VariantContext>();

    private void listFiltersAndExit() {
        out.println("\nAvailable features: " + getAvailableClasses(featureClasses));
        out.println("Available exclusion criteria: " + getAvailableClasses(exclusionClasses) + "\n");
        System.exit(0);
    }

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        featureClasses = PackageUtils.getClassesImplementingInterface(IndependentVariantFeature.class);
        exclusionClasses = PackageUtils.getClassesImplementingInterface(VariantExclusionCriterion.class);

        if (LIST) { listFiltersAndExit(); }

        try {
            sampleNames.put(sampleName.toUpperCase(), "variant");
            vcfHeader = VariantsToVCF.getHeader(this.getToolkit().getArguments(), sampleNames.keySet());
            vcfWriter = new VCFWriter(vcfHeader, VCF_OUT);
            includedWriter = new PrintWriter(INCLUDED_OUT);
            includedWriter.println(GeliTextWriter.headerLine);
            if ( ANNOTATED_OUT != null ) {
                annotatedWriter = new PrintWriter(ANNOTATED_OUT);
                annotatedWriter.print("Chr\tPosition\t");
            }
            requestedFeatures = new ArrayList<IndependentVariantFeature>();
            requestedExclusions = new ArrayList<VariantExclusionCriterion>();

            // Initialize requested features
            if (FEATURES != null) {
                for (String requestedFeatureString : FEATURES) {
                    String[] requestedFeaturePieces = requestedFeatureString.split(":");
                    String requestedFeatureName = requestedFeaturePieces[0];
                    String requestedFeatureArgs = (requestedFeaturePieces.length == 2) ? requestedFeaturePieces[1] : "";

                    boolean found = false;
                    for ( Class featureClass : featureClasses ) {
                        String featureClassName = rationalizeClassName(featureClass);

                        if (requestedFeatureName.equalsIgnoreCase(featureClassName)) {
                            found = true;

                            try {
                                IndependentVariantFeature ivf = (IndependentVariantFeature) featureClass.newInstance();
                                ivf.initialize(requestedFeatureArgs);
                                requestedFeatures.add(ivf);

                                if ( annotatedWriter != null )
                                    annotatedWriter.print(ivf.getStudyHeader() + "\t");
                            } catch (InstantiationException e) {
                                throw new StingException(String.format("Cannot instantiate feature class '%s': must be concrete class", featureClass.getSimpleName()));
                            } catch (IllegalAccessException e) {
                                throw new StingException(String.format("Cannot instantiate feature class '%s': must have no-arg constructor", featureClass.getSimpleName()));
                            }
                        }
                    }

                    if (!found) {
                        throw new StingException("Unknown feature '" + requestedFeatureString + "'.  Issue the '-ls' argument to list available features.");
                    }
                }
            }

            if (EXCLUSIONS != null) {
                for (String requestedExclusionString : EXCLUSIONS) {
                    String[] requestedExclusionPieces = requestedExclusionString.split(":");
                    String requestedExclusionName = requestedExclusionPieces[0];

                    boolean found = false;
                    for ( Class exclusionClass : exclusionClasses ) {
                        String exclusionClassName = rationalizeClassName(exclusionClass);

                        if (requestedExclusionName.equalsIgnoreCase(exclusionClassName)) {
                            found = true;

                            try {
                                HashMap<String,String> requestedArgs = new HashMap<String,String>();
                                if ( requestedExclusionPieces.length == 2 ) {
                                    String[] argStrings = requestedExclusionPieces[1].split(",");
                                    for (int i = 0; i < argStrings.length; i++ ) {
                                        String[] arg = argStrings[i].split("=");
                                        if ( arg.length == 2 )
                                            requestedArgs.put(arg[0], arg[1]);
                                    }
                                }
                                VariantExclusionCriterion vec = (VariantExclusionCriterion) exclusionClass.newInstance();
                                vec.initialize(requestedArgs);
                                requestedExclusions.add(vec);

                                if ( annotatedWriter != null )
                                    annotatedWriter.print(vec.getStudyHeader() + "\t");
                            } catch (InstantiationException e) {
                                throw new StingException(String.format("Cannot instantiate exclusion class '%s': must be concrete class", exclusionClass.getSimpleName()));
                            } catch (IllegalAccessException e) {
                                throw new StingException(String.format("Cannot instantiate exclusion class '%s': must have no-arg constructor", exclusionClass.getSimpleName()));
                            }
                        }
                    }

                    if (!found) {
                        throw new StingException("Unknown exclusion criterion '" + requestedExclusionString + "'.  Issue the '-ls' argument to list available exclusion criteria.");
                    }
                }
            }

            if ( annotatedWriter != null )
                annotatedWriter.print("inDbSNP\tinHapMap\tisHet\n");
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file(s) for writing"));
        }
    }

    /**
     * Trim the 'IVF' or 'VEC' off the feature/exclusion name so the user needn't specify that on the command-line.
     *
     * @param featureClass  the feature class whose name we should rationalize
     * @return  the class name, minus 'IVF'
     */
    private String rationalizeClassName(Class featureClass) {
        String featureClassName = featureClass.getSimpleName();
        String newName = featureClassName.replaceFirst("IVF", "");
        newName = newName.replaceFirst("VEC", "");
        return newName;
    }

    /**
     * Returns a comma-separated list of available classes the user may specify at the command-line.
     *
     * @param classes an ArrayList of classes
     * @return String of available classes 
     */
    private <T> String getAvailableClasses(List<Class<? extends T>> classes) {
        String availableString = "";

        for (int classIndex = 0; classIndex < classes.size(); classIndex++) {
            availableString += rationalizeClassName(classes.get(classIndex)) + (classIndex == classes.size() - 1 ? "" : ",");
        }

        return availableString;
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    /**
     * For each site of interest, rescore the genotype likelihoods by applying the specified feature set.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        RodGeliText variant = (RodGeliText) tracker.lookup("variant", null);

        // Ignore places where we don't have a variant or where the reference base is ambiguous.
        if ( variant == null || BaseUtils.simpleBaseToBaseIndex(ref.getBase()) == -1 )
            return 0;

        VariantContext varContext = new VariantContext(tracker, ref, context, variant);

        // if we're still initializing the context, do so
        if ( windowInitializer != null ) {
            windowInitializer.add(varContext);
            if ( windowInitializer.size() == windowSize ) {
                variantContextWindow = new VariantContextWindow(windowInitializer);
                windowInitializer = null;
            }
        } else {
            variantContextWindow.moveWindow(varContext);
            compute();
        }

        return 1;
    }

    private void compute() {
        // get the current context
        VariantContext context = variantContextWindow.getContext();
        if ( context == null )
            return;
        RodGeliText variant = context.getVariant();

        HashMap<String, Double> exclusionResults = new HashMap<String, Double>();

        if (VERBOSE) { out.println("Original:\n" + variant); }

        GenomeLoc loc = context.getAlignmentContext(true).getLocation();
        if ( annotatedWriter != null )
            annotatedWriter.print(loc.getContig() + "\t" + loc.getStart() + "\t");

        // Apply features that modify the likelihoods and LOD scores
        for ( IndependentVariantFeature ivf : requestedFeatures ) {
            ivf.compute(variantContextWindow);

            double[] weights = ivf.getLikelihoods();

            variant.adjustLikelihoods(weights);

            if (VERBOSE) { out.println(rationalizeClassName(ivf.getClass()) + ":\n  " + variant); }

            if ( annotatedWriter != null )
                annotatedWriter.print(ivf.getStudyInfo() + "\t");
        }

        // Apply exclusion tests that score the variant call
        if (VERBOSE) {
            out.print("InclusionProbabilities:[");
        }

        // Use the filters to score the variant
        String filterFailureString = "";
        double jointInclusionProbability = 1.0;
        for ( VariantExclusionCriterion vec : requestedExclusions ) {
            vec.compute(variantContextWindow);

            String exclusionClassName = rationalizeClassName(vec.getClass());

            Double inclusionProbability = vec.inclusionProbability();
            jointInclusionProbability *= inclusionProbability;
            exclusionResults.put(exclusionClassName, inclusionProbability);

            if (inclusionProbability < INCLUSION_THRESHOLD) {
                filterFailureString += vec.getVCFFilterString() + ";";
            }

            if (VERBOSE) {
                out.print(exclusionClassName + "=" + inclusionProbability + ";");
            }

            if ( annotatedWriter != null )
                annotatedWriter.print(vec.getStudyInfo() + "\t");
        }

        // Decide whether we should keep the call or not
        if (jointInclusionProbability >= INCLUSION_THRESHOLD) {
            includedWriter.println(variant);

            if (VERBOSE) { out.println("] JointInclusionProbability:" + jointInclusionProbability + " State:included\n"); }
        } else {
            if (VERBOSE) { out.println("] JointInclusionProbability:" + jointInclusionProbability + " State:excluded\n"); }
        }

        rodDbSNP dbsnp = (rodDbSNP) context.getTracker().lookup("dbSNP", null);
        if ( annotatedWriter != null ) {
            if ( dbsnp == null )
                annotatedWriter.print("false\tfalse\t");
            else
                annotatedWriter.print(dbsnp.isSNP() + "\t" + dbsnp.isHapmap() + "\t");
            annotatedWriter.println(GenotypeUtils.isHet(variant));
        }

        List<VCFGenotypeRecord> gt = new ArrayList<VCFGenotypeRecord>();
        Map<VCFHeader.HEADER_FIELDS,String> map = new HashMap<VCFHeader.HEADER_FIELDS,String>();
        if ( VariantsToVCF.generateVCFRecord(context.getTracker(), context.getReferenceContext(), context.getAlignmentContext(true), vcfHeader, gt, map, sampleNames, out, false, false) ) {
            if ( !filterFailureString.equals("") )
                map.put(VCFHeader.HEADER_FIELDS.FILTER, filterFailureString);           
            vcfWriter.addRecord(new VCFRecord(vcfHeader, map, "GT:GQ:DP", gt));
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
        // move the window over so that we can filter the last few variants
        if ( windowInitializer != null ) {
            while ( windowInitializer.size() < windowSize )
                windowInitializer.add(null);
            variantContextWindow = new VariantContextWindow(windowInitializer);
        }
        for (int i=0; i < windowSize; i++) {
            variantContextWindow.moveWindow(null);
            compute();
        }

        out.printf("Processed %d loci.\n", result);

        vcfWriter.close();
        includedWriter.close();
        if ( annotatedWriter != null )
            annotatedWriter.close();
    }
}
