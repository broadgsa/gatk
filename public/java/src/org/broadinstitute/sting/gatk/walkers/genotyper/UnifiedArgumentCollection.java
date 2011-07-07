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

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Input;

import java.io.File;


public class UnifiedArgumentCollection {

    // control the various models to be used
    @Argument(fullName = "genotype_likelihoods_model", shortName = "glm", doc = "Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together", required = false)
    public GenotypeLikelihoodsCalculationModel.Model GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;

    @Argument(fullName = "p_nonref_model", shortName = "pnrm", doc = "Non-reference probability calculation model to employ -- EXACT is the default option, while GRID_SEARCH is also available.", required = false)
    public AlleleFrequencyCalculationModel.Model AFmodel = AlleleFrequencyCalculationModel.Model.EXACT;

    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = DiploidSNPGenotypePriors.HUMAN_HETEROZYGOSITY;

    @Argument(fullName = "pcr_error_rate", shortName = "pcr_error", doc = "The PCR error rate to be used for computing fragment-based likelihoods", required = false)
    public Double PCR_error = DiploidSNPGenotypeLikelihoods.DEFAULT_PCR_ERROR_RATE;

    @Argument(fullName = "genotyping_mode", shortName = "gt_mode", doc = "Should we output confident genotypes (i.e. including ref calls) or just the variants?", required = false)
    public GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY;

    @Argument(fullName = "output_mode", shortName = "out_mode", doc = "Should we output confident genotypes (i.e. including ref calls) or just the variants?", required = false)
    public UnifiedGenotyperEngine.OUTPUT_MODE OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY;

    @Argument(fullName = "standard_min_confidence_threshold_for_calling", shortName = "stand_call_conf", doc = "The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be called", required = false)
    public double STANDARD_CONFIDENCE_FOR_CALLING = 30.0;

    @Argument(fullName = "standard_min_confidence_threshold_for_emitting", shortName = "stand_emit_conf", doc = "The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be emitted (and filtered if less than the calling threshold)", required = false)
    public double STANDARD_CONFIDENCE_FOR_EMITTING = 30.0;

    @Argument(fullName = "noSLOD", shortName = "nsl", doc = "If provided, we will not calculate the SLOD", required = false)
    public boolean NO_SLOD = false;


    // control the error modes
    @Hidden
    @Argument(fullName = "assume_single_sample_reads", shortName = "single_sample", doc = "The single sample that we should assume is represented in the input bam (and therefore associate with all reads regardless of whether they have read groups)", required = false)
    public String ASSUME_SINGLE_SAMPLE = null;

    // TODO -- delete me
    @Hidden
    @Argument(fullName = "abort_at_too_much_coverage", doc = "Don't call a site if the downsampled coverage is greater than this value", required = false)
    public int COVERAGE_AT_WHICH_TO_ABORT = -1;


    // control the various parameters to be used
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public int MIN_BASE_QUALTY_SCORE = 17;

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for calling", required = false)
    public int MIN_MAPPING_QUALTY_SCORE = 20;

    @Argument(fullName = "max_deletion_fraction", shortName = "deletions", doc = "Maximum fraction of reads with deletions spanning this locus for it to be callable [to disable, set to < 0 or > 1; default:0.05]", required = false)
    public Double MAX_DELETION_FRACTION = 0.05;


    // indel-related arguments
    @Argument(fullName = "min_indel_count_for_genotyping", shortName = "minIndelCnt", doc = "Minimum number of consensus indels required to trigger genotyping run", required = false)
    public int MIN_INDEL_COUNT_FOR_GENOTYPING = 5;

    @Argument(fullName = "indel_heterozygosity", shortName = "indelHeterozygosity", doc = "Heterozygosity for indel calling", required = false)
    public double INDEL_HETEROZYGOSITY = 1.0/8000;

    @Hidden
    @Argument(fullName = "indelGapContinuationPenalty", shortName = "indelGCP", doc = "Indel gap continuation penalty", required = false)
    public double INDEL_GAP_CONTINUATION_PENALTY = 10.0;

    @Hidden
    @Argument(fullName = "indelGapOpenPenalty", shortName = "indelGOP", doc = "Indel gap open penalty", required = false)
    public double INDEL_GAP_OPEN_PENALTY = 45.0;

    @Hidden
    @Argument(fullName = "indelHaplotypeSize", shortName = "indelHSize", doc = "Indel haplotype size", required = false)
    public int INDEL_HAPLOTYPE_SIZE = 80;
    @Hidden
    @Argument(fullName = "doContextDependentGapPenalties", shortName = "doCDP", doc = "Vary gap penalties by context", required = false)
     public boolean DO_CONTEXT_DEPENDENT_PENALTIES = true;
    //gdebug+
    // experimental arguments, NOT TO BE USED BY ANYONE WHOSE INITIALS AREN'T GDA!!!
    @Hidden
    @Argument(fullName = "getGapPenaltiesFromData", shortName = "dataGP", doc = "Vary gap penalties by context - EXPERIMENTAL, DO NO USE", required = false)
    public boolean GET_GAP_PENALTIES_FROM_DATA = false;

    @Hidden
    @Argument(fullName="indel_recal_file", shortName="recalFile", required=false, doc="Filename for the input covariates table recalibration .csv file - EXPERIMENTAL, DO NO USE")
    public File INDEL_RECAL_FILE = new File("indel.recal_data.csv");

    @Hidden
    @Argument(fullName = "indelDebug", shortName = "indelDebug", doc = "Output indel debug info", required = false)
    public boolean OUTPUT_DEBUG_INDEL_INFO = false;
    @Hidden
    @Argument(fullName = "dovit", shortName = "dovit", doc = "Output indel debug info", required = false)
    public boolean dovit = false;
    @Hidden
    @Argument(fullName = "GSA_PRODUCTION_ONLY", shortName = "GSA_PRODUCTION_ONLY", doc = "don't ever use me", required = false)
    public boolean GSA_PRODUCTION_ONLY = false;
    @Hidden
 
    @Argument(fullName = "exactCalculation", shortName = "exactCalculation", doc = "expt", required = false)
    public ExactAFCalculationModel.ExactCalculation EXACT_CALCULATION_TYPE = ExactAFCalculationModel.ExactCalculation.LINEAR_EXPERIMENTAL;

    @Hidden
     @Argument(fullName = "ignoreSNPAlleles", shortName = "ignoreSNPAlleles", doc = "expt", required = false)
    public boolean IGNORE_SNP_ALLELES = false;


    @Deprecated
    @Argument(fullName="output_all_callable_bases", shortName="all_bases", doc="Please use --output_mode EMIT_ALL_SITES instead" ,required=false)
    private Boolean ALL_BASES_DEPRECATED = false;   

    @Deprecated
    @Argument(fullName="genotype", shortName="genotype", doc="Please use --output_mode EMIT_ALL_CONFIDENT_SITES instead" ,required=false)
    private Boolean GENOTYPE_DEPRECATED = false;


    // Developers must remember to add any newly added arguments to the list here as well otherwise they won't get changed from their default value!
    public UnifiedArgumentCollection clone() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();

        uac.GLmodel = GLmodel;
        uac.AFmodel = AFmodel;
        uac.EXACT_CALCULATION_TYPE = EXACT_CALCULATION_TYPE;
        uac.heterozygosity = heterozygosity;
        uac.PCR_error = PCR_error;
        uac.GenotypingMode = GenotypingMode;
        uac.OutputMode = OutputMode;
        uac.NO_SLOD = NO_SLOD;
        uac.ASSUME_SINGLE_SAMPLE = ASSUME_SINGLE_SAMPLE;
        uac.STANDARD_CONFIDENCE_FOR_CALLING = STANDARD_CONFIDENCE_FOR_CALLING;
        uac.STANDARD_CONFIDENCE_FOR_EMITTING = STANDARD_CONFIDENCE_FOR_EMITTING;
        uac.MIN_BASE_QUALTY_SCORE = MIN_BASE_QUALTY_SCORE;
        uac.MIN_MAPPING_QUALTY_SCORE = MIN_MAPPING_QUALTY_SCORE;
        uac.MAX_DELETION_FRACTION = MAX_DELETION_FRACTION;
        uac.MIN_INDEL_COUNT_FOR_GENOTYPING = MIN_INDEL_COUNT_FOR_GENOTYPING;
        uac.INDEL_HETEROZYGOSITY = INDEL_HETEROZYGOSITY;
        uac.INDEL_GAP_OPEN_PENALTY = INDEL_GAP_OPEN_PENALTY;
        uac.INDEL_GAP_CONTINUATION_PENALTY = INDEL_GAP_CONTINUATION_PENALTY;
        uac.OUTPUT_DEBUG_INDEL_INFO = OUTPUT_DEBUG_INDEL_INFO;
        uac.INDEL_HAPLOTYPE_SIZE = INDEL_HAPLOTYPE_SIZE;
        uac.DO_CONTEXT_DEPENDENT_PENALTIES = DO_CONTEXT_DEPENDENT_PENALTIES;

        uac.GET_GAP_PENALTIES_FROM_DATA = GET_GAP_PENALTIES_FROM_DATA;
        uac.INDEL_RECAL_FILE = INDEL_RECAL_FILE;
        // todo- arguments to remove
        uac.COVERAGE_AT_WHICH_TO_ABORT = COVERAGE_AT_WHICH_TO_ABORT;
        uac.dovit = dovit;
        uac.GSA_PRODUCTION_ONLY = GSA_PRODUCTION_ONLY;
        uac.IGNORE_SNP_ALLELES = IGNORE_SNP_ALLELES;
        
        return uac;
    }


}
