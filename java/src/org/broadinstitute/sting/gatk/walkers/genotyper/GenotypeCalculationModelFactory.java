package org.broadinstitute.sting.gatk.walkers.genotyper;

import static org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel.Model.*;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;

import java.util.Set;


public class GenotypeCalculationModelFactory {
    //private GenotypeCalculationModelFactory() {} // cannot be instantiated

    public static GenotypeCalculationModel.Model getGenotypeCalculationModel(final String name) {
        GenotypeCalculationModel.Model m = valueOf(name);
        if ( m == null )
            throw new RuntimeException("Unexpected GenotypeCalculationModel " + name);
        else
            return m;
    }

    /**
     * General and correct way to create GenotypeCalculationModel objects for arbitrary models
     *
     * @param UAC           The unified argument  collection
     * @param samples       samples in bam
     * @param out           output
     * @return model
     */
    public static GenotypeCalculationModel makeGenotypeCalculation(UnifiedArgumentCollection UAC,
                                                                   Set<String> samples,
                                                                   GenotypeWriter out) {
        switch ( UAC.genotypeModel ) {
            case SINGLE_SAMPLE: return new SingleSampleGenotypeCalculationModel(UAC.baseModel, samples, UAC.defaultPlatform, out, UAC.GENOTYPE, UAC.LOD_THRESHOLD, UAC.MAX_DELETIONS, UAC.VERBOSE);
            //case MULTI_SAMPLE_EM: return new MultiSampleEMGenotypeCalculationModel(UAC.baseModel, samples, UAC.defaultPlatform, out, UAC.GENOTYPE, UAC.LOD_THRESHOLD, UAC.MAX_DELETIONS, UAC.VERBOSE);
            case MULTI_SAMPLE_ALL_MAFS: return new MultiSampleAllMAFsGenotypeCalculationModel(UAC.baseModel, samples, UAC.defaultPlatform, out, UAC.GENOTYPE, UAC.LOD_THRESHOLD, UAC.MAX_DELETIONS, UAC.VERBOSE);
            default: throw new RuntimeException("Unexpected GenotypeCalculationModel " + UAC.genotypeModel);
        }
    }
}