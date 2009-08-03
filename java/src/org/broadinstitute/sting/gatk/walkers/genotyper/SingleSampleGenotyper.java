package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.playground.utils.AlleleMetrics;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.calls.SSGGenotypeCall;

import java.io.File;

@ReadFilters(ZeroMappingQualityReadFilter.class)
public class SingleSampleGenotyper extends LocusWalker<SSGGenotypeCall, GenotypeWriter> {
    // Control output settings
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = true)
    public File VARIANTS_FILE;
    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File to which metrics should be written", required = false)
    public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName = "lod_threshold", shortName = "lod", doc = "The lod threshold on which variants should be filtered", required = false)
    public Double LOD_THRESHOLD = Double.MIN_VALUE;

    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confidient genotypes or just the variants?", required = false)
    public boolean GENOTYPE = false;

    @Argument(fullName = "3BaseErrors", shortName = "3BaseErrors", doc = "Should we use a 3-base error mode (so that P(b_true != b_called | e) == e / 3?", required = false)
    public boolean THREE_BASE_ERRORS = false;

    public enum Caller {
        OLD_AND_BUSTED,
        NEW_HOTNESS
    }

    @Argument(fullName = "caller", doc = "", required = false)
    public Caller caller = Caller.OLD_AND_BUSTED;

    // Control how the genotype hypotheses are weighed
    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false) public Double heterozygosity = GenotypeLikelihoods.HUMAN_HETEROZYGOSITY;
    @Argument(fullName = "priors_hapmap", shortName = "phapmap", doc = "Comma-separated prior likelihoods for Hapmap loci (homref,het,homvar)", required = false) public String PRIORS_HAPMAP = "0.999,1e-3,1e-5";
    @Argument(fullName = "priors_dbsnp", shortName = "pdbsnp", doc = "Comma-separated prior likelihoods for dbSNP loci (homref,het,homvar)", required = false) public String PRIORS_DBSNP = "0.999,1e-3,1e-5";

    @Argument(fullName = "keepQ0Bases", shortName = "keepQ0Bases", doc = "If true, then Q0 bases will be included in the genotyping calculation, and treated as Q1 -- this is really not a good idea", required = false)
    public boolean keepQ0Bases = false;

    public double[] plocus;
    public double[] phapmap;
    public double[] pdbsnp;

    public long nFilteredQ0Bases = 0;

    /**
     * Specify that this walker requires reads.
     *
     * @return true
     */
    public boolean requiresReads() {
        return true;
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
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return (BaseUtils.simpleBaseToBaseIndex(ref) != -1 && context.getReads().size() != 0);
    }

    /**
     * Convert an array string (value1,value2,...,valueN) to a double array
     *
     * @param priorsString the array string of priors
     *
     * @return the same array, but each value is now a double
     */
    private double[] priorsArray(String priorsString) {
        String[] pstrs = priorsString.split(",");
        double[] pdbls = new double[pstrs.length];

        for (int i = 0; i < pstrs.length; i++) {
            pdbls[i] = Double.valueOf(pstrs[i]);
        }

        return pdbls;
    }

    /** Initialize the walker with some sensible defaults */
    public void initialize() {
        plocus = GenotypeLikelihoods.computePriors(heterozygosity);
        phapmap = priorsArray(PRIORS_HAPMAP);
        pdbsnp = priorsArray(PRIORS_DBSNP);
    }

    /**
     * Compute the AlleleFrequencyEstimate at a given locus.
     *
     * @param tracker the meta data tracker
     * @param ref     the reference base
     * @param context contextual information around the locus
     *
     * @return an AlleleFrequencyEstimate object
     */
    public SSGGenotypeCall map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        SSGGenotypeCall oldAndBusted = mapOldAndBusted(tracker, ref, context);
        SSGGenotypeCall newHotness = mapNewHotness(tracker, ref, context);

        if ( ! oldAndBusted.equals(newHotness) ) {
            System.out.printf("Calls not equal:%nold: %s%nnew: %s%n", oldAndBusted, newHotness);
        }

        return caller == Caller.OLD_AND_BUSTED ? oldAndBusted : newHotness;
    }

    private SSGGenotypeCall mapNewHotness(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        GenotypeLikelihoods G = makeGenotypeLikelihood(tracker);
        G.setDiscovery(GENOTYPE); // set it to discovery mode or variant detection mode
        SSGGenotypeCall geno = (SSGGenotypeCall)G.callGenotypes(tracker, ref, pileup);
        return geno;
    }

    private SSGGenotypeCall mapOldAndBusted(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        GenotypeLikelihoods G = makeGenotypeLikelihood(tracker);
        G.setDiscovery(GENOTYPE); // set it to discovery mode or variant detection mode
        SSGGenotypeCall geno = (SSGGenotypeCall)G.callGenotypes(tracker, ref, pileup);
        return geno;
    }    

    /**
     * Calls the underlying, single locus genotype of the sample
     *
     * @param tracker the meta data tracker
     *
     * @return the likelihoods per genotype
     */
    private GenotypeLikelihoods makeGenotypeLikelihood(RefMetaDataTracker tracker) {
        GenotypeLikelihoods G = null;

        if (isHapmapSite(tracker)) {
            G = new GenotypeLikelihoods(GenotypeLikelihoods.HUMAN_HETEROZYGOSITY);
        } else if (isDbSNPSite(tracker)) {
            G = new GenotypeLikelihoods(pdbsnp[0], pdbsnp[1], pdbsnp[2]);
        } else {
            G = new GenotypeLikelihoods(plocus[0], plocus[1], plocus[2]);
        }

        G.filterQ0Bases(! keepQ0Bases);

        return G;
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
    public GenotypeWriter reduceInit() {
        return GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getSAMFileHeader(), VARIANTS_FILE);

    }

    /**
     * If we've found a LOD >= 5 variant, output it to disk.
     *
     * @param call an GenotypeCall object for the variant.
     * @param sum  accumulator for the reduce.
     *
     * @return an empty string
     */
    public GenotypeWriter reduce(SSGGenotypeCall call, GenotypeWriter sum) {
        if (call != null && (GENOTYPE || call.isVariation())) {
            if (call.getConfidenceScore().getScore() > LOD_THRESHOLD)
                sum.addGenotypeCall(call);
        }
        return sum;
    }

    /** Close the variant file. */
    public void onTraversalDone(GenotypeWriter sum) {
        logger.info(String.format("SingleSampleGenotyper filtered %d Q0 bases", nFilteredQ0Bases));
        sum.close();
    }
}


