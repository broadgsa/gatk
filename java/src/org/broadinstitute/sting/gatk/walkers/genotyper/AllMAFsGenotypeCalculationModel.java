package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.List;

public class AllMAFsGenotypeCalculationModel extends GenotypeCalculationModel {

    protected AllMAFsGenotypeCalculationModel() {}

    public List<GenotypeCall> calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // This is just a stub...
        // TODO --- implement me

        // In this alternative approach we calculate the likelihoods of all possible
        // MAFs (1-N) and choose the most likely one.  We can probably be smart and
        // avoid calculating MAFs that we know can't be the most likely...

        return null;
    }
}