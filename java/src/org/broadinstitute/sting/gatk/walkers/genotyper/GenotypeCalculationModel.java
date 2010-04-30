package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;

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
    protected GenotypeWriterFactory.GENOTYPE_FORMAT OUTPUT_FORMAT;
    protected PrintStream verboseWriter;
    protected PrintStream beagleWriter;

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
     * @param outputFormat  output format
     * @param verboseWriter verbose writer
     * @param beagleWriter  beagle writer
     */
    protected void initialize(Set<String> samples,
                              Logger logger,
                              UnifiedArgumentCollection UAC,
                              GenotypeWriterFactory.GENOTYPE_FORMAT outputFormat,
                              PrintStream verboseWriter,
                              PrintStream beagleWriter) {
        this.UAC = UAC.clone();
        this.samples = new TreeSet<String>(samples);
        OUTPUT_FORMAT = outputFormat;
        this.logger = logger;
        this.verboseWriter = verboseWriter;
        this.beagleWriter = beagleWriter;
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
                                                 char ref,
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
                                                         char[] ref,
                                                         GenomeLoc loc,
                                                         Map<String, StratifiedAlignmentContext> stratifiedContexts);

    /**
     * @param tracker   rod data
     *
     * @return the dbsnp rod if there is one at this position
     */
    public static rodDbSNP getDbSNP(RefMetaDataTracker tracker) {
        return rodDbSNP.getFirstRealSNP(tracker.getReferenceMetaData(rodDbSNP.STANDARD_DBSNP_TRACK_NAME));
    }
}