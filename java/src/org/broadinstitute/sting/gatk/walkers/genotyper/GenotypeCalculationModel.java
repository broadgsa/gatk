package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
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
    protected int maxDeletionsInPileup;
    protected String assumedSingleSample;
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
        maxDeletionsInPileup = UAC.MAX_DELETIONS;
        assumedSingleSample = UAC.ASSUME_SINGLE_SAMPLE;
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
     * @param tracker   rod data
     * @param ref       reference base
     * @param context   alignment context
     * @param priors    priors to use for GL
     *
     * @return calls
     */
    public abstract Pair<List<Genotype>, GenotypeLocusData> calculateGenotype(RefMetaDataTracker tracker,
                                                                             char ref,
                                                                             AlignmentContext context,
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

    /**
     * Create the mapping from sample to alignment contexts.
     * @param context   original alignment context
     *
     * @return stratified contexts, or null if there are too many deletions at this position
     */
    protected HashMap<String, AlignmentContextBySample> splitContextBySample(AlignmentContext context) {

        HashMap<String, AlignmentContextBySample> contexts = new HashMap<String, AlignmentContextBySample>();

        int deletionsInPileup = 0;
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        for (int i = 0; i < reads.size(); i++) {

            // get the read and offset
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            // find the sample
            String sample;
            SAMReadGroupRecord readGroup = read.getReadGroup();
            if ( readGroup == null ) {
                if ( assumedSingleSample == null )
                    throw new StingException("Missing read group for read " + read.getReadName());
                sample = assumedSingleSample;
            } else {
                sample = readGroup.getSample();
            }

            // create a new context object if this is the first time we're seeing a read for this sample
            AlignmentContextBySample myContext = contexts.get(sample);
            if ( myContext == null ) {
                myContext = new AlignmentContextBySample(context.getLocation());
                contexts.put(sample, myContext);
            }

            // check for deletions
            if ( offset == -1 ) {
                // are there too many deletions in the pileup?
                if ( ++deletionsInPileup > maxDeletionsInPileup && maxDeletionsInPileup >= 0 )
                    return null;
            }

            // add the read to this sample's context
            // note that bad bases are added to the context (for DoC calculations later)
            myContext.add(read, offset);
        }

        return contexts;
    }

    
    /**
     * A class to keep track of the alignment context observed for a given sample.
     * we currently store the overall context and strand-stratified sets,
     * but any arbitrary stratification could be added.
    */
    protected enum StratifiedContext { OVERALL, FORWARD, REVERSE }
    protected class AlignmentContextBySample {

        private AlignmentContext overall = null;
        private AlignmentContext forward = null;
        private AlignmentContext reverse = null;
        private GenomeLoc loc;

        private ArrayList<SAMRecord> allReads = new ArrayList<SAMRecord>();
        private ArrayList<SAMRecord> forwardReads = new ArrayList<SAMRecord>();
        private ArrayList<SAMRecord> reverseReads = new ArrayList<SAMRecord>();

        private ArrayList<Integer> allOffsets = new ArrayList<Integer>();
        private ArrayList<Integer> forwardOffsets = new ArrayList<Integer>();
        private ArrayList<Integer> reverseOffsets = new ArrayList<Integer>();


        AlignmentContextBySample(GenomeLoc loc) {
            this.loc = loc;
        }

        public AlignmentContext getContext(StratifiedContext context) {
            switch ( context ) {
                case OVERALL: return getOverallContext();
                case FORWARD: return getForwardContext();
                case REVERSE: return getReverseContext();
            }
            return null;
        }

        private AlignmentContext getOverallContext() {
            if ( overall == null )
                overall = new AlignmentContext(loc, allReads, allOffsets);
            return overall;
        }

        private AlignmentContext getForwardContext() {
            if ( forward == null )
                forward = new AlignmentContext(loc, forwardReads, forwardOffsets);
            return forward;
        }

        private AlignmentContext getReverseContext() {
            if ( reverse == null )
                reverse = new AlignmentContext(loc, reverseReads, reverseOffsets);
            return reverse;
        }

        public void add(SAMRecord read, int offset) {
            if ( read.getReadNegativeStrandFlag() ) {
                reverseReads.add(read);
                reverseOffsets.add(offset);
            } else {
                forwardReads.add(read);
                forwardOffsets.add(offset);
            }
            allReads.add(read);
            allOffsets.add(offset);
         }
    }
}