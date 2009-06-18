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
    @Argument(fullName="features", shortName="F", doc="Comma-separated list of feature tests to apply to genotype posteriors.") public String FEATURES;
    @Argument(fullName="variants_out", shortName="VO", doc="File to which modified variants should be written") public File VARIANTS_OUT;

    private String[] features;
    private PrintWriter vwriter;

    public void initialize() {
        features = FEATURES.split(",");

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

        for (String feature : features) {
            IndependentVariantFeature ivf;
            
            if (feature.equalsIgnoreCase("binomial")) { ivf = new IVFBinomialStrand(); }
            else { throw new StingException(String.format("Cannot understand feature '%s'\n", feature)); }

            variant.multiplyLikelihoods(ivf.compute(ref, context));
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
