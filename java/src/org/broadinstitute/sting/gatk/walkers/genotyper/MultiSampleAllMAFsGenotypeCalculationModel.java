package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.Set;

public class MultiSampleAllMAFsGenotypeCalculationModel extends GenotypeCalculationModel {

    public MultiSampleAllMAFsGenotypeCalculationModel(BaseMismatchModel baseModel,
                                                      Set<String> samples,
                                                      EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform platform,
                                                      GenotypeWriter out,
                                                      boolean genotypeMode,
                                                      double lod,
                                                      int maxDeletions,
                                                      boolean verbose) {
        super(baseModel, samples, platform, out, genotypeMode, lod, maxDeletions, verbose);
    }


    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public boolean calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // This is just a stub...
        // TODO --- implement me

        // In this alternative approach we calculate the likelihoods of all possible
        // MAFs (1-N) and choose the most likely one.  We can probably be smart and
        // avoid calculating MAFs that we know can't be the most likely...

        return true;
    }
}