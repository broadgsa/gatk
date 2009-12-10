package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

/**
 * The model representing how we calculate a genotype given the priors and a pile
 * of bases and quality scores
 */
public abstract class GenotypeCalculationModel implements Cloneable {

    public enum Model {
        EM_POINT_ESTIMATE,
        JOINT_ESTIMATE,
        POOLED
    }

    protected BaseMismatchModel baseModel;
    protected Set<String> samples;
    protected Logger logger;
    protected double heterozygosity;
    protected EmpiricalSubstitutionProbabilities.SequencerPlatform defaultPlatform;
    protected GenotypeWriterFactory.GENOTYPE_FORMAT OUTPUT_FORMAT;
    protected boolean ALL_BASE_MODE;
    protected boolean GENOTYPE_MODE;
    protected int POOL_SIZE;
    protected double CONFIDENCE_THRESHOLD;
    protected double MINIMUM_ALLELE_FREQUENCY;
    protected boolean REPORT_SLOD;
    protected PrintWriter verboseWriter;

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
        this.samples = new TreeSet<String>(samples);
        this.logger = logger;
        baseModel = UAC.baseModel;
        heterozygosity = UAC.heterozygosity;
        defaultPlatform = UAC.defaultPlatform;
        OUTPUT_FORMAT = UAC.VAR_FORMAT;
        ALL_BASE_MODE = UAC.ALL_BASES;
        GENOTYPE_MODE = UAC.GENOTYPE;
        POOL_SIZE = UAC.POOLSIZE;
        CONFIDENCE_THRESHOLD = UAC.CONFIDENCE_THRESHOLD;
        MINIMUM_ALLELE_FREQUENCY = UAC.MINIMUM_ALLELE_FREQUENCY;
        if ( UAC.VERBOSE != null ) {
            try {
                verboseWriter = new PrintWriter(UAC.VERBOSE);
                initializeVerboseWriter(verboseWriter);
            } catch (FileNotFoundException e) {
                throw new StingException("Could not open file " + UAC.VERBOSE + " for writing");
            }
        }
        REPORT_SLOD = ! UAC.NO_SLOD;
    }

    protected void initializeVerboseWriter(PrintWriter writer) { };

    public void close() {
        if ( verboseWriter != null )
            verboseWriter.close();
    }

    /**
     * Must be overridden by concrete subclasses
     * @param tracker              rod data
     * @param ref                  reference base
     * @param loc                  GenomeLoc
     * @param stratifiedContexts   stratified alignment contexts
     * @param priors               priors to use for GL
     *
     * @return calls
     */
    public abstract Pair<VariationCall, List<Genotype>> calculateGenotype(RefMetaDataTracker tracker,
                                                                          char ref,
                                                                          GenomeLoc loc,
                                                                          Map<String, StratifiedAlignmentContext> stratifiedContexts,
                                                                          DiploidGenotypePriors priors);

    /**
     * @param tracker   rod data
     *
     * @return the dbsnp rod if there is one at this position
     */
    public static rodDbSNP getDbSNP(RefMetaDataTracker tracker) {
        return rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null));
    }

    /**
     * Determine whether we're at a Hapmap site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a Hapmap site, false if otherwise
     */
    private static boolean isHapmapSite(RefMetaDataTracker tracker) {
        return tracker.getTrackData("hapmap", null) != null;
    }
}