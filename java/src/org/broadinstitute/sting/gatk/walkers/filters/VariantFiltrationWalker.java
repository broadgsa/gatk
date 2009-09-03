package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.geli.GeliTextWriter;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;


/**
 * VariantFiltrationWalker applies specified conditionally independent features and filters to pre-called variants.
 * The former modifiesthe likelihoods of each genotype, while the latter makes a decision to include or exclude a
 * variant outright.  At the moment, the variants are expected to be in gelitext format.
 */
@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=rodVariants.class))
public class VariantFiltrationWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="variants_out_head", shortName="VOH", doc="File to which modified variants should be written") public String VARIANTS_OUT_HEAD;
    @Argument(fullName="features", shortName="F", doc="Feature test (optionally with arguments) to apply to genotype posteriors.  Syntax: 'testname[:arguments]'", required=false) public String[] FEATURES;
    @Argument(fullName="exclusion_criterion", shortName="X", doc="Exclusion test (optionally with arguments) to apply to variant call.  Syntax: 'testname[:arguments]'", required=false) public String[] EXCLUSIONS;
    @Argument(fullName="inclusion_threshold", shortName="IT", doc="The product of the probability to include variants based on these filters must be greater than the value specified here in order to be included", required=false) public Double INCLUSION_THRESHOLD = 0.9;
    @Argument(fullName="verbose", shortName="V", doc="Show how the variant likelihoods are changing with the application of each feature") public Boolean VERBOSE = false;
    @Argument(fullName="list", shortName="ls", doc="List the available features and exclusion criteria and exit") public Boolean LIST = false;

    private List<Class<? extends IndependentVariantFeature>> featureClasses;
    private List<Class<? extends VariantExclusionCriterion>> exclusionClasses;

    private PrintWriter variantsWriter;
    private PrintWriter paramsWriter;
    private HashMap<String, PrintWriter> exclusionWriters;

    private ArrayList<IndependentVariantFeature> requestedFeatures;
    private ArrayList<VariantExclusionCriterion> requestedExclusions;

    // the structures necessary to initialize and maintain a windowed context
    private VariantContextWindow variantContextWindow;
    private static final int windowSize = 10;  // 10 variants on either end of the current one
    private ArrayList<VariantContext> windowInitializer = new ArrayList<VariantContext>();

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        featureClasses = PackageUtils.getClassesImplementingInterface(IndependentVariantFeature.class);
        exclusionClasses = PackageUtils.getClassesImplementingInterface(VariantExclusionCriterion.class);

        if (LIST) {
            out.println("\nAvailable features: " + getAvailableClasses(featureClasses));
            out.println("Available exclusion criteria: " + getAvailableClasses(exclusionClasses) + "\n");
            System.exit(0);
        }

        try {
            variantsWriter = new PrintWriter(VARIANTS_OUT_HEAD + ".included.geli.calls");
            variantsWriter.println(GeliTextWriter.headerLine);

            paramsWriter = new PrintWriter(VARIANTS_OUT_HEAD + ".params.out");
            paramsWriter.print("Chr\tPosition\t");

            requestedFeatures = new ArrayList<IndependentVariantFeature>();
            requestedExclusions = new ArrayList<VariantExclusionCriterion>();

            // Initialize requested features
            if (FEATURES != null) {
                for (String requestedFeatureString : FEATURES) {
                    String[] requestedFeaturePieces = requestedFeatureString.split(":");
                    String requestedFeatureName = requestedFeaturePieces[0];
                    String requestedFeatureArgs = (requestedFeaturePieces.length == 2) ? requestedFeaturePieces[1] : "";

                    for ( Class featureClass : featureClasses ) {
                        String featureClassName = rationalizeClassName(featureClass);

                        if (requestedFeatureName.equalsIgnoreCase(featureClassName)) {
                            try {
                                IndependentVariantFeature ivf = (IndependentVariantFeature) featureClass.newInstance();
                                ivf.initialize(requestedFeatureArgs);
                                requestedFeatures.add(ivf);

                                paramsWriter.print(ivf.getStudyHeader() + "\t");
                            } catch (InstantiationException e) {
                                throw new StingException(String.format("Cannot instantiate feature class '%s': must be concrete class", featureClass.getSimpleName()));
                            } catch (IllegalAccessException e) {
                                throw new StingException(String.format("Cannot instantiate feature class '%s': must have no-arg constructor", featureClass.getSimpleName()));
                            }
                        }
                    }
                }
            }

            // Initialize requested exclusion criteria
            exclusionWriters = new HashMap<String, PrintWriter>();

            if (EXCLUSIONS != null) {
                for (String requestedExclusionString : EXCLUSIONS) {
                    String[] requestedExclusionPieces = requestedExclusionString.split(":");
                    String requestedExclusionName = requestedExclusionPieces[0];
                    String requestedExclusionArgs = (requestedExclusionPieces.length == 2) ? requestedExclusionPieces[1] : "";

                    for ( Class exclusionClass : exclusionClasses ) {
                        String exclusionClassName = rationalizeClassName(exclusionClass);

                        if (requestedExclusionName.equalsIgnoreCase(exclusionClassName)) {
                            try {
                                VariantExclusionCriterion vec = (VariantExclusionCriterion) exclusionClass.newInstance();
                                vec.initialize(requestedExclusionArgs);
                                requestedExclusions.add(vec);

                                paramsWriter.print(vec.getStudyHeader() + "\t");

                                PrintWriter writer = new PrintWriter(VARIANTS_OUT_HEAD + ".excluded." + exclusionClassName + ".geli.calls");
                                writer.println(GeliTextWriter.headerLine);

                                exclusionWriters.put(exclusionClassName, writer);
                            } catch (InstantiationException e) {
                                throw new StingException(String.format("Cannot instantiate exclusion class '%s': must be concrete class", exclusionClass.getSimpleName()));
                            } catch (IllegalAccessException e) {
                                throw new StingException(String.format("Cannot instantiate exclusion class '%s': must have no-arg constructor", exclusionClass.getSimpleName()));
                            }
                        }
                    }
                }
            }

            paramsWriter.print("inDbSNP\tinHapMap\tisHet\n");
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
        rodVariants variant = (rodVariants) tracker.lookup("variant", null);

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
        rodVariants variant = context.getVariant();

        HashMap<String, Double> exclusionResults = new HashMap<String, Double>();

        if (VERBOSE) { out.println("Original:\n" + variant); }

        GenomeLoc loc = context.getAlignmentContext(true).getLocation();
        paramsWriter.print(loc.getContig() + "\t" + loc.getStart() + "\t");

        // Apply features that modify the likelihoods and LOD scores
        for ( IndependentVariantFeature ivf : requestedFeatures ) {
            ivf.compute(variantContextWindow);

            double[] weights = ivf.getLikelihoods();

            variant.adjustLikelihoods(weights);

            if (VERBOSE) { out.println(rationalizeClassName(ivf.getClass()) + ":\n  " + variant); }

            paramsWriter.print(ivf.getStudyInfo() + "\t");
        }

        // Apply exclusion tests that score the variant call
        if (VERBOSE) {
            out.print("InclusionProbabilities:[");
        }

        // Use the filters to score the variant
        double jointInclusionProbability = 1.0;
        for ( VariantExclusionCriterion vec : requestedExclusions ) {
            vec.compute(variantContextWindow);

            String exclusionClassName = rationalizeClassName(vec.getClass());

            Double inclusionProbability = vec.inclusionProbability();
            jointInclusionProbability *= inclusionProbability;
            exclusionResults.put(exclusionClassName, inclusionProbability);

            if (inclusionProbability < INCLUSION_THRESHOLD) {
                PrintWriter ewriter = exclusionWriters.get(exclusionClassName);
                if (ewriter != null) {
                    ewriter.println(variant);
                    ewriter.flush();
                }
            }

            if (VERBOSE) {
                out.print(exclusionClassName + "=" + inclusionProbability + ";");
            }

            paramsWriter.print(vec.getStudyInfo() + "\t");
        }

        // Decide whether we should keep the call or not
        if (jointInclusionProbability >= INCLUSION_THRESHOLD) {
            variantsWriter.println(variant);

            if (VERBOSE) { out.println("] JointInclusionProbability:" + jointInclusionProbability + " State:included\n"); }
        } else {
            if (VERBOSE) { out.println("] JointInclusionProbability:" + jointInclusionProbability + " State:excluded\n"); }
        }

        rodDbSNP dbsnp = (rodDbSNP) context.getTracker().lookup("dbSNP", null);
        if ( dbsnp == null ) {
            paramsWriter.print("false\tfalse\t");
        } else {
            paramsWriter.print(dbsnp.isSNP() + "\t" + dbsnp.isHapmap() + "\t");
        }

        paramsWriter.println(GenotypeUtils.isHet(variant));
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

        variantsWriter.close();
        paramsWriter.close();

        for (PrintWriter ewriter : exclusionWriters.values()) {
            ewriter.close();
        }
    }
}
