package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.GenotypeMetaData;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.apache.log4j.Logger;

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


    /**
     * Create the mapping from sample to alignment contexts.
     * @param context   original alignment context
     *
     * @return stratified contexts, or null if there are too many deletions at this position
     */
    protected HashMap<String, AlignmentContextBySample> splitContextBySample(AlignmentContext context) {

        HashMap<String, AlignmentContextBySample> contexts = new HashMap<String, AlignmentContextBySample>();

        int deletionsInPile = 0;
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        for (int i = 0; i < reads.size(); i++) {

            // get the read and offset
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            // find the sample; special case for pools
            String sample;
            if ( POOLED_INPUT ) {
                sample = "POOL";
            } else {
                SAMReadGroupRecord readGroup = read.getReadGroup();
                if ( readGroup == null ) {
                    if ( assumedSingleSample == null )
                        throw new StingException("Missing read group for read " + read.getReadName());
                    sample = assumedSingleSample;
                } else {
                    sample = readGroup.getSample();
                }
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
                if ( ++deletionsInPile > maxDeletionsInPileup )
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