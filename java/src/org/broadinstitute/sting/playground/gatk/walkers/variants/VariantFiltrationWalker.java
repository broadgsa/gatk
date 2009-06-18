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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * VariantFiltrationWalker applies specified conditionally independent features to pre-called variants, thus modifying
 * the likelihoods of each genotype.  At the moment, the variants are expected to be in gelitext format.
 */
@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=rodVariants.class))
public class VariantFiltrationWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="features", shortName="F", doc="Feature test (optionally with arguments) to apply to genotype posteriors.  Syntax: 'testname[:arguments]'") public String[] FEATURES;
    @Argument(fullName="variants_out", shortName="VO", doc="File to which modified variants should be written") public File VARIANTS_OUT;
    @Argument(fullName="verbose", shortName="V", doc="Show how the variant likelihoods are changing with the application of each feature") public Boolean VERBOSE = false;

    private ArrayList<Class> featureClasses;
    private PrintWriter vwriter;

    /**
     * Trim the 'IVF' off the feature name so the user needn't specify that on the command-line.
     *
     * @param featureClass  the feature class whose name we should rationalize
     * @return  the class name, minus 'IVF'
     */
    private String rationalizeFeatureClassName(Class featureClass) {
        String featureClassName = featureClass.getSimpleName();
        return featureClassName.replaceFirst("IVF", "");
    }

    /**
     * Returns a comma-separated list of available features the user may specify at the command-line.
     *
     * @return String of available features
     */
    private String getAvailableFeatureClasses() {
        String featureString = "";

        for (int featureClassIndex = 0; featureClassIndex < featureClasses.size(); featureClassIndex++) {
            featureString += rationalizeFeatureClassName(featureClasses.get(featureClassIndex)) + (featureClassIndex == featureClasses.size() - 1 ? "" : ",");
        }

        return featureString;
    }

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        try {
            vwriter = new PrintWriter(VARIANTS_OUT);
            vwriter.println(AlleleFrequencyEstimate.geliHeaderString());

            featureClasses = PackageUtils.getClassesImplementingInterface(IndependentVariantFeature.class);
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file '%s' for writing", VARIANTS_OUT.getAbsolutePath()));
        }
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

            for (String requestedFeatureString : FEATURES) {
                String[] requestedFeaturePieces = requestedFeatureString.split(":");
                String requestedFeatureName = requestedFeaturePieces[0];
                String requestedFeatureArgs = (requestedFeaturePieces.length == 2) ? requestedFeaturePieces[1] : "";

                int notYetSeenFeature = 0;
                for ( Class featureClass : featureClasses ) {
                    String featureClassName = rationalizeFeatureClassName(featureClass);

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
                    throw new StingException(String.format("Unknown feature '%s'. Valid features are '%s'", requestedFeatureName, getAvailableFeatureClasses()));
                }

                if (VERBOSE) { System.out.println(); }
            }

            vwriter.println(variant);

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
    }
}
