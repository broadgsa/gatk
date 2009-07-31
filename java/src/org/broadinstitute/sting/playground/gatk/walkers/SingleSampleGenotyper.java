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
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = true) public File VARIANTS_FILE;
    @Argument(fullName = "metrics_out", shortName = "metout", doc = "File to which metrics should be written", required = false) public File METRICS_FILE = new File("/dev/stderr");
    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File to which metrics should be written", required = false) public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName = "lod_threshold", shortName = "lod", doc = "The lod threshold on which variants should be filtered", required = false)public Double LOD_THRESHOLD = Double.MIN_VALUE;
    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confidient genotypes or just the variants?", required = false) public boolean GENOTYPE = false;

    @Argument(fullName = "3BaseErrors", shortName = "3BaseErrors", doc = "Should we use a 3-base error mode (so that P(b_true != b_called | e) == e / 3?", required = false) public boolean THREE_BASE_ERRORS = false;

    // Control periodic reporting features
    @Argument(fullName = "metrics_interval", shortName = "metint", doc = "Number of loci to process between metrics reports", required = false) public Integer METRICS_INTERVAL = 50000;
    @Argument(fullName = "suppress_metrics", shortName = "printmets", doc = "If specified, don't display metrics", required = false) public Boolean SUPPRESS_METRICS = false;

    // Control what features we use in calling variants
    //@Argument(fullName = "ignore_secondary_bases", shortName = "nosb", doc = "Ignore secondary base examination", required = false) public Boolean IGNORE_SECONDARY_BASES = false;
    @Argument(fullName = "call_indels", shortName = "indels", doc = "Call indels as well as point mutations", required = false) public Boolean CALL_INDELS = false;

    // Control how the genotype hypotheses are weighed
    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false) public Double heterozygosity = 0.001;
    @Argument(fullName = "priors_hapmap", shortName = "phapmap", doc = "Comma-separated prior likelihoods for Hapmap loci (homref,het,homvar)", required = false) public String PRIORS_HAPMAP = "0.999,1e-3,1e-5";
    @Argument(fullName = "priors_dbsnp", shortName = "pdbsnp", doc = "Comma-separated prior likelihoods for dbSNP loci (homref,het,homvar)", required = false) public String PRIORS_DBSNP = "0.999,1e-3,1e-5";
    @Argument(fullName = "priors_2nd_on", shortName = "p2ndon", doc = "Comma-separated prior likelihoods for the secondary bases of primary on-genotype bases (AA,AC,AG,AT,CC,CG,CT,GG,GT,TT)", required = false) public String PRIORS_2ND_ON = "0.000,0.302,0.366,0.142,0.000,0.548,0.370,0.000,0.319,0.000";
    @Argument(fullName = "priors_2nd_off", shortName = "p2ndoff", doc = "Comma-separated prior likelihoods for the secondary bases of primary off-genotype bases (AA,AC,AG,AT,CC,CG,CT,GG,GT,TT)", required = false) public String PRIORS_2ND_OFF = "0.480,0.769,0.744,0.538,0.575,0.727,0.768,0.589,0.762,0.505";

    // Control various sample-level settings
    @Argument(fullName = "sample_name_regex", shortName = "sample_name_regex", doc = "Replaces the sample name specified in the BAM read group with the value supplied here", required = false) public String SAMPLE_NAME_REGEX = null;

    @Argument(fullName = "keepQ0Bases", shortName = "keepQ0Bases", doc = "If true, then Q0 bases will be included in the genotyping calculation, and treated as Q1 -- this is really not a good idea", required = false)
    public boolean keepQ0Bases = false;

    public AlleleMetrics metricsOut;
    public String sampleName;

    public double[] plocus;
    public double[] phapmap;
    public double[] pdbsnp;
    public double[] p2ndon;
    public double[] p2ndoff;

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

    private double[] computePriors(double h) {
        double[] pdbls = new double[3];
        pdbls[0] = 1.0 - (3.0 * h / 2.0);
        pdbls[1] = h;
        pdbls[2] = h / 2.0;
        return pdbls;
    }

    /** Initialize the walker with some sensible defaults */
    public void initialize() {
        metricsOut = new AlleleMetrics(METRICS_FILE, LOD_THRESHOLD);

        plocus = computePriors(heterozygosity);
        phapmap = priorsArray(PRIORS_HAPMAP);
        pdbsnp = priorsArray(PRIORS_DBSNP);
        p2ndon = priorsArray(PRIORS_2ND_ON);
        p2ndoff = priorsArray(PRIORS_2ND_OFF);
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
        rationalizeSampleName(context.getReads().get(0));
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        GenotypeLikelihoods G = callGenotype(tracker);
        G.setDiscovery(GENOTYPE); // set it to discovery mode or variant detection mode
        SSGGenotypeCall geno = (SSGGenotypeCall)G.callGenotypes(tracker, ref, pileup);
        if (geno != null) {
            metricsOut.nextPosition(geno, tracker);
        }
        if (!SUPPRESS_METRICS) {
            metricsOut.printMetricsAtLocusIntervals(METRICS_INTERVAL);
        }
        return geno;
    }

    /**
     * Sometimes the sample names in the BAM files get screwed up.  Fix it here if we can.
     *
     * @param read a read from the pileup (assuming all the reads have the same sample name)
     *
     * @return a repaired sample name
     */
    private String rationalizeSampleName(SAMRecord read) {
        String RG = (String) (read.getAttribute("RG"));
        SAMReadGroupRecord read_group_record = read.getHeader().getReadGroup(RG);

        if (read_group_record != null) {
            String localSampleName = read.getHeader().getReadGroup(RG).getSample();
            if (SAMPLE_NAME_REGEX != null) {
                localSampleName = localSampleName.replaceAll(SAMPLE_NAME_REGEX, "$1");
            }
            if (sampleName == null) {
                sampleName = localSampleName;
            } else {
                if (!sampleName.equals(localSampleName)) {
                    throw new StingException(String.format("Samples appear to have been mixed up: expected '%s' but found '%s'.", sampleName, localSampleName));
                }
            }
        }

        return sampleName;
    }

    /**
     * Calls the underlying, single locus genotype of the sample
     *
     * @param tracker the meta data tracker
     *
     * @return the likelihoods per genotype
     */
    private GenotypeLikelihoods callGenotype(RefMetaDataTracker tracker) {
        GenotypeLikelihoods G = null;

        if (isHapmapSite(tracker)) {
            G = new GenotypeLikelihoods(THREE_BASE_ERRORS, phapmap[0], phapmap[1], phapmap[2], p2ndon, p2ndoff, keepQ0Bases);
        } else if (isDbSNPSite(tracker)) {
            G = new GenotypeLikelihoods(THREE_BASE_ERRORS, pdbsnp[0], pdbsnp[1], pdbsnp[2], p2ndon, p2ndoff, keepQ0Bases);
        } else {
            G = new GenotypeLikelihoods(THREE_BASE_ERRORS, plocus[0], plocus[1], plocus[2], p2ndon, p2ndoff, keepQ0Bases);
        }
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


