/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardCallerArgumentCollection;
import org.broadinstitute.sting.utils.pairhmm.PairHMM;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

public class UnifiedArgumentCollection extends StandardCallerArgumentCollection {

    @Argument(fullName = "genotype_likelihoods_model", shortName = "glm", doc = "Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together", required = false)
    public GenotypeLikelihoodsCalculationModel.Model GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;

    /**
     * The PCR error rate is independent of the sequencing error rate, which is necessary because we cannot necessarily
     * distinguish between PCR errors vs. sequencing errors.  The practical implication for this value is that it
     * effectively acts as a cap on the base qualities.
     */
    @Argument(fullName = "pcr_error_rate", shortName = "pcr_error", doc = "The PCR error rate to be used for computing fragment-based likelihoods", required = false)
    public Double PCR_error = DiploidSNPGenotypeLikelihoods.DEFAULT_PCR_ERROR_RATE;

    /**
     * Note that calculating the SLOD increases the runtime by an appreciable amount.
     */
    @Argument(fullName = "computeSLOD", shortName = "slod", doc = "If provided, we will calculate the SLOD (SB annotation)", required = false)
    public boolean COMPUTE_SLOD = false;

    /**
     * Depending on the value of the --max_alternate_alleles argument, we may genotype only a fraction of the alleles being sent on for genotyping.
     * Using this argument instructs the genotyper to annotate (in the INFO field) the number of alternate alleles that were originally discovered at the site.
     */
    @Argument(fullName = "annotateNDA", shortName = "nda", doc = "If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site", required = false)
    public boolean ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = false;

    /**
     * The PairHMM implementation to use for -glm INDEL genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    @Argument(fullName = "pair_hmm_implementation", shortName = "pairHMM", doc = "The PairHMM implementation to use for -glm INDEL genotype likelihood calculations", required = false)
    public PairHMM.HMM_IMPLEMENTATION pairHMM = PairHMM.HMM_IMPLEMENTATION.ORIGINAL;

    /**
     * The minimum confidence needed in a given base for it to be used in variant calling.  Note that the base quality of a base
     * is capped by the mapping quality so that bases on reads with low mapping quality may get filtered out depending on this value.
     * Note too that this argument is ignored in indel calling. In indel calling, low-quality ends of reads are clipped off (with fixed threshold of Q20).
     */
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public int MIN_BASE_QUALTY_SCORE = 17;

    @Argument(fullName = "max_deletion_fraction", shortName = "deletions", doc = "Maximum fraction of reads with deletions spanning this locus for it to be callable [to disable, set to < 0 or > 1; default:0.05]", required = false)
    public Double MAX_DELETION_FRACTION = 0.05;

    // indel-related arguments
    /**
     * A candidate indel is genotyped (and potentially called) if there are this number of reads with a consensus indel at a site.
     * Decreasing this value will increase sensitivity but at the cost of larger calling time and a larger number of false positives.
     */
    @Argument(fullName = "min_indel_count_for_genotyping", shortName = "minIndelCnt", doc = "Minimum number of consensus indels required to trigger genotyping run", required = false)
    public int MIN_INDEL_COUNT_FOR_GENOTYPING = 5;

    /**
     * Complementary argument to minIndelCnt.  Only samples with at least this fraction of indel-containing reads will contribute
     * to counting and overcoming the threshold minIndelCnt.  This parameter ensures that in deep data you don't end
     * up summing lots of super rare errors up to overcome the 5 read default threshold.  Should work equally well for
     * low-coverage and high-coverage samples, as low coverage samples with any indel containing reads should easily over
     * come this threshold.
     */
    @Argument(fullName = "min_indel_fraction_per_sample", shortName = "minIndelFrac", doc = "Minimum fraction of all reads at a locus that must contain an indel (of any allele) for that sample to contribute to the indel count for alleles", required = false)
    public double MIN_INDEL_FRACTION_PER_SAMPLE = 0.25;

    /**
     * This argument informs the prior probability of having an indel at a site.
     */
    @Argument(fullName = "indel_heterozygosity", shortName = "indelHeterozygosity", doc = "Heterozygosity for indel calling", required = false)
    public double INDEL_HETEROZYGOSITY = 1.0/8000;

    @Advanced
    @Argument(fullName = "indelGapContinuationPenalty", shortName = "indelGCP", doc = "Indel gap continuation penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10", required = false)
    public byte INDEL_GAP_CONTINUATION_PENALTY = 10;

    @Advanced
    @Argument(fullName = "indelGapOpenPenalty", shortName = "indelGOP", doc = "Indel gap open penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10", required = false)
    public byte INDEL_GAP_OPEN_PENALTY = 45;

    @Hidden
    @Argument(fullName = "indelHaplotypeSize", shortName = "indelHSize", doc = "Indel haplotype size", required = false)
    public int INDEL_HAPLOTYPE_SIZE = 80;

    @Hidden
    @Argument(fullName = "indelDebug", shortName = "indelDebug", doc = "Output indel debug info", required = false)
    public boolean OUTPUT_DEBUG_INDEL_INFO = false;

    @Hidden
    @Argument(fullName = "ignoreSNPAlleles", shortName = "ignoreSNPAlleles", doc = "expt", required = false)
    public boolean IGNORE_SNP_ALLELES = false;

    /*
        Generalized ploidy argument (debug only): squash all reads into a single pileup without considering sample info
     */
    @Hidden
    @Argument(fullName = "allReadsSP", shortName = "dl", doc = "expt", required = false)
    public boolean TREAT_ALL_READS_AS_SINGLE_POOL = false;

    /*
       Generalized ploidy argument (debug only): When building site error models, ignore lane information and build only
       sample-level error model
     */
    @Argument(fullName = "ignoreLaneInfo", shortName = "ignoreLane", doc = "Ignore lane when building error model, error model is then per-site", required = false)
    public boolean IGNORE_LANE_INFO = false;

    /*
        Generalized ploidy argument: VCF file that contains truth calls for reference sample. If a reference sample is included through argument -refsample,
        then this argument is required.
     */
    @Input(fullName="reference_sample_calls", shortName = "referenceCalls", doc="VCF file with the truth callset for the reference sample", required=false)
    RodBinding<VariantContext> referenceSampleRod;

    /*
        Reference sample name: if included, a site-specific error model will be built in order to improve calling quality. This requires ideally
        that a bar-coded reference sample be included with the polyploid/pooled data in a sequencing experimental design.
        If argument is absent, no per-site error model is included and calling is done with a generalization of traditional statistical calling.
     */
    @Argument(shortName="refsample", fullName="reference_sample_name", doc="Reference sample name.", required=false)
    String referenceSampleName;

    /*
        Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
     */
    @Argument(shortName="ploidy", fullName="sample_ploidy", doc="Plody (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).", required=false)
    public int samplePloidy = VariantContextUtils.DEFAULT_PLOIDY;

    @Hidden
    @Argument(shortName="minqs", fullName="min_quality_score", doc="Min quality score to consider. Smaller numbers process faster. Default: Q1.", required=false)
    byte minQualityScore= 1;

    @Hidden
    @Argument(shortName="maxqs", fullName="max_quality_score", doc="Max quality score to consider. Smaller numbers process faster. Default: Q40.", required=false)
    byte maxQualityScore= 40;

    @Hidden
    @Argument(shortName="site_prior", fullName="site_quality_prior", doc="Phred-Scaled prior quality of the site. Default: Q20.", required=false)
    byte phredScaledPrior = 20;

    @Hidden
    @Argument(shortName = "min_call_power", fullName = "min_power_threshold_for_calling", doc="The minimum confidence in the error model to make a call. Number should be between 0 (no power requirement) and 1 (maximum power required).", required = false)
    double minPower = 0.95;

    @Hidden
    @Argument(shortName = "min_depth", fullName = "min_reference_depth", doc="The minimum depth required in the reference sample in order to make a call.", required = false)
    int minReferenceDepth = 100;

    @Hidden
    @Argument(shortName="ef", fullName="exclude_filtered_reference_sites", doc="Don't include in the analysis sites where the reference sample VCF is filtered. Default: false.", required=false)
    boolean EXCLUDE_FILTERED_REFERENCE_SITES = false;

    /**
     * Create a new UAC with defaults for all UAC arguments
     */
    public UnifiedArgumentCollection() {
        super();
    }

    /**
     * Create a new UAC based on the information only our in super-class scac and defaults for all UAC arguments
     * @param scac
     */
    public UnifiedArgumentCollection(final StandardCallerArgumentCollection scac) {
        super(scac);
    }

    /**
     * Create a new UAC with all parameters having the values in uac
     *
     * @param uac
     */
    public UnifiedArgumentCollection(final UnifiedArgumentCollection uac) {
        // Developers must remember to add any newly added arguments to the list here as well otherwise they won't get changed from their default value!
        super(uac);

        this.GLmodel = uac.GLmodel;
        this.AFmodel = uac.AFmodel;
        this.PCR_error = uac.PCR_error;
        this.COMPUTE_SLOD = uac.COMPUTE_SLOD;
        this.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = uac.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED;
        this.MIN_BASE_QUALTY_SCORE = uac.MIN_BASE_QUALTY_SCORE;
        this.MAX_DELETION_FRACTION = uac.MAX_DELETION_FRACTION;
        this.MIN_INDEL_COUNT_FOR_GENOTYPING = uac.MIN_INDEL_COUNT_FOR_GENOTYPING;
        this.MIN_INDEL_FRACTION_PER_SAMPLE = uac.MIN_INDEL_FRACTION_PER_SAMPLE;
        this.INDEL_HETEROZYGOSITY = uac.INDEL_HETEROZYGOSITY;
        this.INDEL_GAP_OPEN_PENALTY = uac.INDEL_GAP_OPEN_PENALTY;
        this.INDEL_GAP_CONTINUATION_PENALTY = uac.INDEL_GAP_CONTINUATION_PENALTY;
        this.OUTPUT_DEBUG_INDEL_INFO = uac.OUTPUT_DEBUG_INDEL_INFO;
        this.INDEL_HAPLOTYPE_SIZE = uac.INDEL_HAPLOTYPE_SIZE;
        this.TREAT_ALL_READS_AS_SINGLE_POOL = uac.TREAT_ALL_READS_AS_SINGLE_POOL;
        this.referenceSampleRod = uac.referenceSampleRod;
        this.referenceSampleName = uac.referenceSampleName;
        this.samplePloidy = uac.samplePloidy;
        this.maxQualityScore = uac.minQualityScore;
        this.phredScaledPrior = uac.phredScaledPrior;
        this.minPower = uac.minPower;
        this.minReferenceDepth = uac.minReferenceDepth;
        this.EXCLUDE_FILTERED_REFERENCE_SITES = uac.EXCLUDE_FILTERED_REFERENCE_SITES;
        this.IGNORE_LANE_INFO = uac.IGNORE_LANE_INFO;
        this.pairHMM = uac.pairHMM;

        // todo- arguments to remove
        this.IGNORE_SNP_ALLELES = uac.IGNORE_SNP_ALLELES;
    }
}
