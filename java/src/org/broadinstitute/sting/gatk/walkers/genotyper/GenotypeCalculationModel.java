package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.PrintStream;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;


/**
 * The model representing how we calculate a genotype given the priors and a pile
 * of bases and quality scores
 */
public abstract class GenotypeCalculationModel implements Cloneable {

    public enum Model {
        JOINT_ESTIMATE,
        JOINT_ESTIMATE_EXPT_GL,
        INDELS
    }

    protected UnifiedArgumentCollection UAC;
    protected Set<String> samples;
    protected Logger logger;
    protected PrintStream verboseWriter;

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
     * @param verboseWriter verbose writer
     */
    protected void initialize(Set<String> samples,
                              Logger logger,
                              UnifiedArgumentCollection UAC,
                              PrintStream verboseWriter) {
        this.UAC = UAC.clone();
        this.samples = new TreeSet<String>(samples);
        this.logger = logger;
        this.verboseWriter = verboseWriter;
    }

    /**
     * Must be overridden by concrete subclasses
     * @param tracker              rod data
     * @param ref                  reference base
     * @param loc                  GenomeLoc
     * @param stratifiedContexts   stratified alignment contexts
     * @param priors               priors to use for GL
     *
     * @return call
     */
    public abstract VariantCallContext callLocus(RefMetaDataTracker tracker,
                                                 byte ref,
                                                 GenomeLoc loc,
                                                 Map<String, StratifiedAlignmentContext> stratifiedContexts,
                                                 DiploidGenotypePriors priors);

    /**
     * @param tracker              rod data
     * @param ref                  reference base
     * @param loc                  GenomeLoc
     * @param stratifiedContexts   stratified alignment contexts
     *
     * @return call
     */
    public abstract VariantCallContext callExtendedLocus(RefMetaDataTracker tracker,
                                                         byte[] ref,
                                                         GenomeLoc loc,
                                                         Map<String, StratifiedAlignmentContext> stratifiedContexts);

    /**
     * @param tracker   rod data
     *
     * @return the dbsnp rod if there is one at this position
     */
    public static DbSNPFeature getDbSNP(RefMetaDataTracker tracker) {
        return DbSNPHelper.getFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
    }
}