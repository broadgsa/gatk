package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.playground.utils.AlleleMetrics;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import org.broadinstitute.sting.playground.utils.IndelLikelihood;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import org.broadinstitute.sting.utils.genotype.glf.GLFRecord;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

@ReadFilters(ZeroMappingQualityReadFilter.class)
public class SingleSampleGenotyper extends LocusWalker<AlleleFrequencyEstimate, String> {
    // Control output settings
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = true) public File VARIANTS_FILE;
    @Argument(fullName = "metrics_out", shortName = "metout", doc = "File to which metrics should be written", required = false) public File METRICS_FILE = new File("/dev/stderr");
    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File to which metrics should be written", required = false) public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName = "lod_threshold", shortName = "lod", doc = "The lod threshold on which variants should be filtered", required = false)public Double LOD_THRESHOLD = Double.MIN_VALUE;
    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confidient genotypes or just the variants?", required = false)
    public boolean GENOTYPE = false;

    // Control periodic reporting features
    @Argument(fullName = "metrics_interval", shortName = "metint", doc = "Number of loci to process between metrics reports", required = false) public Integer METRICS_INTERVAL = 50000;
    @Argument(fullName = "suppress_metrics", shortName = "printmets", doc = "If specified, don't display metrics", required = false) public Boolean SUPPRESS_METRICS = false;

    // Control what features we use in calling variants
    @Argument(fullName = "ignore_secondary_bases", shortName = "nosb", doc = "Ignore secondary base examination", required = false) public Boolean IGNORE_SECONDARY_BASES = false;
    @Argument(fullName = "call_indels", shortName = "indels", doc = "Call indels as well as point mutations", required = false) public Boolean CALL_INDELS = false;

    // Control how the genotype hypotheses are weighed
    @Argument(fullName = "priors_any_locus", shortName = "plocus", doc = "Comma-separated prior likelihoods for any locus (homref,het,homvar)", required = false) public String PRIORS_ANY_LOCUS = "0.999,1e-3,1e-5";
    @Argument(fullName = "priors_hapmap", shortName = "phapmap", doc = "Comma-separated prior likelihoods for Hapmap loci (homref,het,homvar)", required = false) public String PRIORS_HAPMAP = "0.999,1e-3,1e-5";
    @Argument(fullName = "priors_dbsnp", shortName = "pdbsnp", doc = "Comma-separated prior likelihoods for dbSNP loci (homref,het,homvar)", required = false) public String PRIORS_DBSNP = "0.999,1e-3,1e-5";
    @Argument(fullName = "priors_2nd_on", shortName = "p2ndon", doc = "Comma-separated prior likelihoods for the secondary bases of primary on-genotype bases (AA,AC,AG,AT,CC,CG,CT,GG,GT,TT)", required = false) public String PRIORS_2ND_ON = "0.000,0.302,0.366,0.142,0.000,0.548,0.370,0.000,0.319,0.000";
    @Argument(fullName = "priors_2nd_off", shortName = "p2ndoff", doc = "Comma-separated prior likelihoods for the secondary bases of primary off-genotype bases (AA,AC,AG,AT,CC,CG,CT,GG,GT,TT)", required = false) public String PRIORS_2ND_OFF = "0.480,0.769,0.744,0.538,0.575,0.727,0.768,0.589,0.762,0.505";

    // Control various sample-level settings
    @Argument(fullName = "sample_name_regex", shortName = "sample_name_regex", doc = "Replaces the sample name specified in the BAM read group with the value supplied here", required = false) public String SAMPLE_NAME_REGEX = null;

    public AlleleMetrics metricsOut;
    public PrintStream variantsOut;
    public String sampleName;
    private GenotypeWriter mGenotypeWriter;

    public double[] plocus;
    public double[] phapmap;
    public double[] pdbsnp;
    public double[] p2ndon;
    public double[] p2ndoff;

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
        metricsOut = new AlleleMetrics(METRICS_FILE, LOD_THRESHOLD);
        if (this.VAR_FORMAT == GenotypeWriterFactory.GENOTYPE_FORMAT.GLF) {
            mGenotypeWriter = GenotypeWriterFactory.create(GenotypeWriterFactory.GENOTYPE_FORMAT.GLF, GenomeAnalysisEngine.instance.getEngine().getSAMHeader(), VARIANTS_FILE);
        } else {
            try {
                variantsOut = new PrintStream(VARIANTS_FILE);
		    } catch (FileNotFoundException e) {
                err.format("Unable to open file '%s'. Perhaps the parent directory does not exist or is read-only.\n", VARIANTS_FILE.getAbsolutePath());
			    System.exit(-1);
		    }
            if (this.VAR_FORMAT == GenotypeWriterFactory.GENOTYPE_FORMAT.GELI) {
                variantsOut.println(AlleleFrequencyEstimate.geliHeaderString());
            } else if (this.VAR_FORMAT == GenotypeWriterFactory.GENOTYPE_FORMAT.TABULAR) {
                variantsOut.println(AlleleFrequencyEstimate.asTabularStringHeader());
            } else {
                throw new StingException("Unsupported single sample genotyper output format: " + this.VAR_FORMAT.toString());
            }
        }
        plocus = priorsArray(PRIORS_ANY_LOCUS);
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
    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        rationalizeSampleName(context.getReads().get(0));

        AlleleFrequencyEstimate freq = getAlleleFrequency(tracker, Character.toUpperCase(ref), context, sampleName);

        if (freq != null) {
            metricsOut.nextPosition(freq, tracker);
        }
        if (!SUPPRESS_METRICS) {
            metricsOut.printMetricsAtLocusIntervals(METRICS_INTERVAL);
        }

        return freq;
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
     * Compute the allele frequency of the underlying genotype at the given locus.
     *
     * @param tracker     the meta data tracker
     * @param ref         the reference base
     * @param context     contextual information around the locus
     * @param sample_name the name of the sample
     *
     * @return the allele frequency estimate
     */
    private AlleleFrequencyEstimate getAlleleFrequency(RefMetaDataTracker tracker, char ref, LocusContext context, String sample_name) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        // Handle indels, but don't do anything with the result yet.
        IndelLikelihood I = (CALL_INDELS) ? callIndel(context, reads, offsets) : null;

        // Handle single-base polymorphisms.
        GenotypeLikelihoods G = callGenotype(tracker, ref, pileup, reads, offsets);

        return G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods, sample_name);
    }

    /**
     * Calls the underlying, single locus genotype of the sample
     *
     * @param tracker the meta data tracker
     * @param ref     the reference base
     * @param pileup  the pileup object for the given locus
     * @param reads   the reads that overlap this locus
     * @param offsets the offsets per read that identify the base at this locus
     *
     * @return the likelihoods per genotype
     */
    private GenotypeLikelihoods callGenotype(RefMetaDataTracker tracker, char ref, ReadBackedPileup pileup, List<SAMRecord> reads, List<Integer> offsets) {
        GenotypeLikelihoods G;

        if (isHapmapSite(tracker)) {
            G = new GenotypeLikelihoods(phapmap[0], phapmap[1], phapmap[2], p2ndon, p2ndoff);
        } else if (isDbSNPSite(tracker)) {
            G = new GenotypeLikelihoods(pdbsnp[0], pdbsnp[1], pdbsnp[2], p2ndon, p2ndoff);
        } else {
            G = new GenotypeLikelihoods(plocus[0], plocus[1], plocus[2], p2ndon, p2ndoff);
        }

        for (int i = 0; i < reads.size(); i++) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
        }

        G.ApplyPrior(ref, this.alt_allele, this.allele_frequency_prior);
        if (!IGNORE_SECONDARY_BASES && pileup.getBases().length() < 750) {
            G.applySecondBaseDistributionPrior(pileup.getBases(), pileup.getSecondaryBasePileup());
        }

        return G;
    }

    /**
     * Compute the likelihood of an indel at this locus.
     *
     * @param context contextual information around the locus
     * @param reads   the reads that overlap this locus
     * @param offsets the offsets per read that identify the base at this locus
     *
     * @return the likelihood of the indel at this location
     */
    private IndelLikelihood callIndel(LocusContext context, List<SAMRecord> reads, List<Integer> offsets) {
        String[] indels = BasicPileup.indelPileup(reads, offsets);
        IndelLikelihood indelCall = new IndelLikelihood(indels, 1e-4);

        if (!indelCall.getType().equals("ref")) {
            System.out.printf("INDEL %s %s\n", context.getLocation(), indelCall);
        }

        return indelCall;
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

    double allele_frequency_prior = -1;
    char alt_allele;

    /**
     * Accessor for PoolCaller to set the allele frequency prior for this sample.
     *
     * @param freq the allele frequency
     * @param alt  the alternate allele
     */
    public void setAlleleFrequencyPrior(double freq, char alt) {
        this.allele_frequency_prior = freq;
        this.alt_allele = alt;
    }

    /**
     * Initialize values appropriately for the reduce step.
     *
     * @return an empty string
     */
    public String reduceInit() {
        return "";
    }

    /**
     * If we've found a LOD >= 5 variant, output it to disk.
     *
     * @param alleleFreq an AlleleFrequencyEstimage object for the variant.
     * @param sum        accumulator for the reduce.
     *
     * @return an empty string
     */
    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) {
        //System.out.printf("AlleleFreqEstimate %s %f %f %f%n", alleleFreq,alleleFreq.lodVsNextBest, alleleFreq.lodVsRef, LOD_THRESHOLD );
        if (alleleFreq != null && (GENOTYPE ? alleleFreq.lodVsNextBest : alleleFreq.lodVsRef) >= LOD_THRESHOLD) {
            if (this.VAR_FORMAT == GenotypeWriterFactory.GENOTYPE_FORMAT.GELI) {
                variantsOut.println(alleleFreq.asGeliString());
            } else if (this.VAR_FORMAT == GenotypeWriterFactory.GENOTYPE_FORMAT.TABULAR) {
                variantsOut.println(alleleFreq.asTabularString());
            } else if (this.VAR_FORMAT == GenotypeWriterFactory.GENOTYPE_FORMAT.GLF) {
                SAMSequenceRecord rec = GenomeLocParser.getContigInfo(alleleFreq.location.getContig());

                double[] likelihoods = new double[10];
                for (int x = 0; x < alleleFreq.posteriors.length; x++) {
                    likelihoods[x] = GLFRecord.LIKELIHOOD_SCALE_FACTOR * alleleFreq.genotypeLikelihoods.likelihoods[x];
                }

                LikelihoodObject obj = new LikelihoodObject(likelihoods, LikelihoodObject.LIKELIHOOD_TYPE.LOG);

                obj.setLikelihoodType(LikelihoodObject.LIKELIHOOD_TYPE.NEGITIVE_LOG);
                this.mGenotypeWriter.addGenotypeCall(rec,(int)alleleFreq.location.getStart(),0.0f,alleleFreq.ref,alleleFreq.depth,obj);
            }
        }

        return "";
    }

    /** Close the variant file. */
    public void onTraversalDone(String sum) { 
        if (this.VAR_FORMAT == GenotypeWriterFactory.GENOTYPE_FORMAT.GLF) {
            mGenotypeWriter.close();
        } else {
            this.variantsOut.close();
        }
    }
}
