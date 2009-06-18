package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;

import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;

@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=rodVariants.class))
public class VariantFiltrationWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="features", shortName="F", doc="Feature test (optionally with arguments) to apply to genotype posteriors.  Syntax: 'testname:arg1,arg2,...,argN'") public String[] FEATURES;
    @Argument(fullName="variants_out", shortName="VO", doc="File to which modified variants should be written") public File VARIANTS_OUT;

    private PrintWriter vwriter;

    public void initialize() {
        try {
            vwriter = new PrintWriter(VARIANTS_OUT);

            vwriter.println(AlleleFrequencyEstimate.geliHeaderString());
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file '%s' for writing", VARIANTS_OUT.getAbsolutePath()));
        }
    }

    public Integer reduceInit() { return 0; }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        rodVariants variant = (rodVariants) tracker.lookup("variant", null);

        for (String feature : FEATURES) {
            String[] featurePieces = feature.split(":");
            String featureName = featurePieces[0];
            String featureArgs = featurePieces[1];

            IndependentVariantFeature ivf;
            
            if (featureName.equalsIgnoreCase("binomialstrand")) { ivf = new IVFBinomialStrand(featureArgs); }
            else { throw new StingException(String.format("Cannot understand feature '%s'", featureName)); }

            variant.adjustLikelihoods(ivf.compute(ref, context));
            vwriter.println(variant);
        }

        return 1;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + 1;
    }

    public void onTraversalDone(Integer result) {
        out.printf("Processed %d loci.\n", result);

        vwriter.close();
    }
}
