package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.PackageUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.gatk.walkers.variants.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * VariantFiltrationWalker applies specified conditionally independent features to pre-called variants, thus modifying
 * the likelihoods of each genotype.  At the moment, the variants are expected to be in gelitext format.
 */
@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=rodVariants.class))
public class VariantFiltrationWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="variants_out_head", shortName="VOH", doc="File to which modified variants should be written") public String VARIANTS_OUT_HEAD;
    @Argument(fullName="features", shortName="F", doc="Feature test (optionally with arguments) to apply to genotype posteriors.  Syntax: 'testname[:arguments]'", required=false) public String[] FEATURES;
    @Argument(fullName="exclusion_criterion", shortName="X", doc="Exclusion test (optionally with arguments) to apply to variant call.  Syntax: 'testname[:arguments]'", required=false) public String[] EXCLUSIONS;
    @Argument(fullName="verbose", shortName="V", doc="Show how the variant likelihoods are changing with the application of each feature") public Boolean VERBOSE = false;
    @Argument(fullName="list", shortName="ls", doc="List the available features and exclusion criteria and exit") public Boolean LIST = false;

    private ArrayList<Class> featureClasses;
    private ArrayList<Class> exclusionClasses;
    private PrintWriter vwriter;
    private HashMap<String, PrintWriter> ewriters;

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
            vwriter = new PrintWriter(VARIANTS_OUT_HEAD + ".included.geli.calls");
            vwriter.println(AlleleFrequencyEstimate.geliHeaderString());

            ewriters = new HashMap<String, PrintWriter>();

            if (EXCLUSIONS != null) {
                for (String requestedExclusionString : EXCLUSIONS) {
                    String[] requestedExclusionPieces = requestedExclusionString.split(":");
                    String requestedExclusionName = requestedExclusionPieces[0];

                    for ( Class exclusionClass : exclusionClasses ) {
                        String exclusionClassName = rationalizeClassName(exclusionClass);

                        if (requestedExclusionName.equalsIgnoreCase(exclusionClassName)) {
                            PrintWriter writer = new PrintWriter(VARIANTS_OUT_HEAD + ".excluded." + exclusionClassName + ".geli.calls");
                            writer.println(AlleleFrequencyEstimate.geliHeaderString());

                            ewriters.put(exclusionClassName, writer);
                        }
                    }
                }
            }
        } catch (FileNotFoundException e) {
            //throw new StingException(String.format("Could not open file '%s' for writing", VARIANTS_OUT.getAbsolutePath()));
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
    private String getAvailableClasses(ArrayList<Class> classes) {
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
    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        rodVariants variant = (rodVariants) tracker.lookup("variant", null);

        // Ignore places where we don't have a variant or where the reference base is ambiguous.
        if (variant != null && BaseUtils.simpleBaseToBaseIndex(ref) != -1) {
            if (VERBOSE) { out.println("Original:\n  " + variant); }

            // Apply features that modify the likelihoods and LOD scores
            if (FEATURES != null) {
                for (String requestedFeatureString : FEATURES) {
                    String[] requestedFeaturePieces = requestedFeatureString.split(":");
                    String requestedFeatureName = requestedFeaturePieces[0];
                    String requestedFeatureArgs = (requestedFeaturePieces.length == 2) ? requestedFeaturePieces[1] : "";

                    int notYetSeenFeature = 0;
                    for ( Class featureClass : featureClasses ) {
                        String featureClassName = rationalizeClassName(featureClass);

                        if (requestedFeatureName.equalsIgnoreCase(featureClassName)) {
                            try {
                                IndependentVariantFeature ivf = (IndependentVariantFeature) featureClass.newInstance();
                                ivf.initialize(requestedFeatureArgs);

                                variant.adjustLikelihoods(ivf.compute(ref, context));

                                if (VERBOSE) { out.println(featureClassName + ":\n  " + variant); }
                            } catch (InstantiationException e) {
                                throw new StingException(String.format("Cannot instantiate feature class '%s': must be concrete class", featureClass.getSimpleName()));
                            } catch (IllegalAccessException e) {
                                throw new StingException(String.format("Cannot instantiate feature class '%s': must have no-arg constructor", featureClass.getSimpleName()));
                            }
                        } else {
                            notYetSeenFeature++;
                        }
                    }

                    if (notYetSeenFeature == featureClasses.size()) {
                        throw new StingException(String.format("Unknown feature '%s'. Valid features are '%s'", requestedFeatureName, getAvailableClasses(featureClasses)));
                    }
                }
            }

            // Apply exclusion tests that accept or reject the variant call
            ArrayList<String> exclusionResults = new ArrayList<String>();

            if (EXCLUSIONS != null) {
                for (String requestedExclusionString : EXCLUSIONS) {
                    String[] requestedExclusionPieces = requestedExclusionString.split(":");
                    String requestedExclusionName = requestedExclusionPieces[0];
                    String requestedExclusionArgs = (requestedExclusionPieces.length == 2) ? requestedExclusionPieces[1] : "";

                    int notYetSeenExclusion = 0;
                    for ( Class exclusionClass : exclusionClasses ) {
                        String exclusionClassName = rationalizeClassName(exclusionClass);

                        if (requestedExclusionName.equalsIgnoreCase(exclusionClassName)) {
                            try {
                                VariantExclusionCriterion vec = (VariantExclusionCriterion) exclusionClass.newInstance();
                                vec.initialize(requestedExclusionArgs);

                                boolean excludeResult = vec.exclude(ref, context, variant);

                                if (excludeResult) {
                                    exclusionResults.add(exclusionClassName);
                                }
                            } catch (InstantiationException e) {
                                throw new StingException(String.format("Cannot instantiate exclusion class '%s': must be concrete class", exclusionClass.getSimpleName()));
                            } catch (IllegalAccessException e) {
                                throw new StingException(String.format("Cannot instantiate exclusion class '%s': must have no-arg constructor", exclusionClass.getSimpleName()));
                            }
                        } else {
                            notYetSeenExclusion++;
                        }
                    }

                    if (notYetSeenExclusion == exclusionClasses.size()) {
                        throw new StingException(String.format("Unknown exclusion '%s'. Valid exclusions are '%s'", requestedExclusionName, getAvailableClasses(exclusionClasses)));
                    }
                }
            }

            if (exclusionResults.size() > 0) {
                String exclusions = "";
                for (int i = 0; i < exclusionResults.size(); i++) {
                    exclusions += exclusionResults.get(i) + (i == exclusionResults.size() - 1 ? "" : ",");

                    PrintWriter writer = ewriters.get(exclusionResults.get(i));
                    if (writer != null) {
                        writer.println(variant);
                    }
                }
                
                if (VERBOSE) {
                    out.printf("Exclusions: %s\n", exclusions);
                }
            } else {
                vwriter.println(variant);
            }

            if (VERBOSE) { out.println(); }

            return 1;
        }

        return 0;
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

        vwriter.close();

        for (PrintWriter writer : ewriters.values()) {
            writer.close();
        }
    }
}
