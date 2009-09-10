package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;

import java.io.File;

@ReadFilters(ZeroMappingQualityReadFilter.class)
public class SingleSampleGenotyper extends LocusWalker<SSGenotypeCall, SingleSampleGenotyper.CallResult> {
    // Control output settings
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = false)
    public File VARIANTS_FILE = null;

    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File format to be used", required = false)
    public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName = "lod_threshold", shortName = "lod", doc = "The lod threshold on which variants should be filtered", required = false)
    public double LOD_THRESHOLD = Double.MIN_VALUE;

    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confident genotypes or just the variants?", required = false)
    public boolean GENOTYPE = false;

    // Control how the genotype hypotheses are weighed
    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = DiploidGenotypePriors.HUMAN_HETEROZYGOSITY;

    @Argument(fullName = "baseModel", shortName = "m", doc = "Base substitution model to employ -- EMPIRICAL is the recommended default, but it's possible to select the ONE_STATE and THREE_STATE models for comparison purposes", required = false)
    public BaseMismatchModel baseModel = BaseMismatchModel.EMPIRICAL;

    @Argument(fullName = "verbose", shortName = "v", doc = "EXPERIMENTAL", required = false)
    public boolean VERBOSE = false;

    @Argument(fullName = "platform", shortName = "pl", doc = "Causes the genotyper to assume that reads without PL header TAG are this platform.  Defaults to null, indicating that the system will throw a runtime exception when such reads are detected", required = false)
    public EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform defaultPlatform = null;

    @Argument(fullName = "disableCache", doc = "[ADVANCED] If true, we won't use the caching system.  This argument is for testing purposes only", required = false)
    public boolean disableCache = false;


    public class CallResult {
        long nConfidentCalls = 0;
        long nNonConfidentCalls = 0;
        long nCalledBases = 0;
        GenotypeWriter writer;

        CallResult(GenotypeWriter writer) {
            this.writer = writer;
        }
        
        public String toString() {
            return String.format("SSG: %d confident and %d non-confident calls were made at %d bases", 
                    nConfidentCalls, nNonConfidentCalls, nCalledBases);
        }
    }

     /**
     * Filter out loci to ignore (at an ambiguous base in the reference or a locus with zero coverage).
     *
     * @param tracker the meta data tracker
     * @param ref     the reference base
     * @param context contextual information around the locus
     *
     * @return true if we should look at this locus, false otherwise
     */
    public boolean filter(RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        return (BaseUtils.simpleBaseToBaseIndex(ref) != -1 && context.getReads().size() != 0);
    }

    /** Initialize the walker with some sensible defaults */
    public void initialize() {
        //GenotypeLikelihoods.clearCache();
        // nothing to do
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param context contextual information around the locus
     */
    public SSGenotypeCall map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext context) {
        char ref = refContext.getBase();
        if ( BaseUtils.isRegularBase(ref) ) {
            DiploidGenotypePriors priors = new DiploidGenotypePriors(ref, heterozygosity, DiploidGenotypePriors.PROB_OF_TRISTATE_GENOTYPE);

            // setup GenotypeLike object
            GenotypeLikelihoods gl = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);

            gl.setVerbose(VERBOSE);
            gl.setEnableCacheFlag(! disableCache);
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
            gl.add(pileup, true);
            gl.validate();

            return new SSGenotypeCall(context.getLocation(), ref,gl, pileup);
        } else {
            return null;
        }
    }

    /**
     * Determine whether we're at a Hapmap site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a Hapmap site, false if otherwise
     */
    private boolean isHapmapSite(RefMetaDataTracker tracker) {
        return tracker.lookup("hapmap", null) != null;
    }

    /**
     * Determine whether we're at a dbSNP site
     *
     * @param tracker the meta data tracker
     *
     * @return true if we're at a dbSNP site, false if otherwise
     */
    private boolean isDbSNPSite(RefMetaDataTracker tracker) {
        return tracker.lookup("dbsnp", null) != null;
    }

    /**
     * Initialize values appropriately for the reduce step.
     *
     * @return an empty string
     */
    public CallResult reduceInit() {
        if ( VARIANTS_FILE != null )
            return new CallResult(GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), VARIANTS_FILE));
        else
            return new CallResult(GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), out));
    }

    /**
     * If we've found a LOD variant or callable genotype, output it to disk.
     *
     * @param call an GenotypeCall object for the variant.
     * @param sum  accumulator for the reduce.
     *
     * @return an empty string
     */
    public CallResult reduce(SSGenotypeCall call, CallResult sum) {
        sum.nCalledBases++;

        if (call != null && (GENOTYPE || call.isVariant(call.getReference()))) {
            double confidence = (GENOTYPE) ? call.getNegLog10PError() : call.toVariation().getNegLog10PError();
            if (confidence >= LOD_THRESHOLD) {
                sum.nConfidentCalls++;
                //System.out.printf("Call %s%n", call);
                sum.writer.addGenotypeCall(call);
            } else
                sum.nNonConfidentCalls++;
        }
        return sum;
    }

    /** Close the variant file. */
    public void onTraversalDone(CallResult sum) {
        sum.writer.close();
    }
}


