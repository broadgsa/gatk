package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.GenotypeMetaData;
import org.broadinstitute.sting.utils.Pair;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * The model representing how we calculate a genotype given the priors and a pile
 * of bases and quality scores
 */
public abstract class GenotypeCalculationModel implements Cloneable {

    public enum Model {
        EM_POINT_ESTIMATE,
        EM_ALL_MAFS
    }

    protected BaseMismatchModel baseModel;
    protected Set<String> samples;
    protected Logger logger;
    protected double heterozygosity;
    protected EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform defaultPlatform;
    protected boolean GENOTYPE_MODE;
    protected boolean POOLED_INPUT;
    protected int POOL_SIZE;
    protected double LOD_THRESHOLD;
    protected int maxDeletionsInPileup;
    protected String assumedSingleSample;
    protected boolean VERBOSE;

    /**
     * Create a new GenotypeCalculationModel object
     */
    protected GenotypeCalculationModel() {}


    /**
     * Initialize the GenotypeCalculationModel object
     * Assumes that out is not null
     * @param samples       samples in input bam
     * @param logger        logger
     * @param UAC           unified arg collection
     */
    protected void initialize(Set<String> samples,
                              Logger logger,
                              UnifiedArgumentCollection UAC) {
        this.samples = samples;
        this.logger = logger;
        baseModel = UAC.baseModel;
        heterozygosity = UAC.heterozygosity;
        defaultPlatform = UAC.defaultPlatform;
        GENOTYPE_MODE = UAC.GENOTYPE;
        POOLED_INPUT = UAC.POOLED;
        POOL_SIZE = UAC.POOLSIZE;
        LOD_THRESHOLD = UAC.LOD_THRESHOLD;
        maxDeletionsInPileup = UAC.MAX_DELETIONS;
        assumedSingleSample = UAC.ASSUME_SINGLE_SAMPLE;
        VERBOSE = UAC.VERBOSE;
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        GenotypeCalculationModel gcm = (GenotypeCalculationModel)super.clone();
        gcm.samples = new HashSet<String>(samples);
        gcm.logger = logger;
        gcm.baseModel = baseModel;
        gcm.heterozygosity = heterozygosity;
        gcm.defaultPlatform = defaultPlatform;
        gcm.GENOTYPE_MODE = GENOTYPE_MODE;
        gcm.POOLED_INPUT = POOLED_INPUT;
        gcm.LOD_THRESHOLD = LOD_THRESHOLD;
        gcm.maxDeletionsInPileup = maxDeletionsInPileup;
        gcm.assumedSingleSample = assumedSingleSample;
        gcm.VERBOSE = VERBOSE;
        return gcm;
    }

    public void setAssumedSingleSample(String sample) {
        assumedSingleSample = sample;
    }

    /**
     * Must be overridden by concrete subclasses
     * @param tracker   rod data
     * @param ref       reference base
     * @param context   alignment context
     * @param priors    priors to use for GL
     *
     * @return calls
     */
    public abstract Pair<List<GenotypeCall>, GenotypeMetaData> calculateGenotype(RefMetaDataTracker tracker,
                                                                                 char ref,
                                                                                 AlignmentContext context,
                                                                                 DiploidGenotypePriors priors);
}