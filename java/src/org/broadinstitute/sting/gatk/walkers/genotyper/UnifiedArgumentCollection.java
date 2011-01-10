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


public class UnifiedArgumentCollection {

    // control the various models to be used
    @Argument(fullName = "genotype_likelihoods_model", shortName = "glm", doc = "Genotype likelihoods calculation model to employ -- SNP is the default option, while DINDEL is also available for calling indels.", required = false)
    public GenotypeLikelihoodsCalculationModel.Model GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;

    @Argument(fullName = "p_nonref_model", shortName = "pnrm", doc = "Non-reference probability calculation model to employ -- EXACT is the default option, while GRID_SEARCH is also available.", required = false)
    public AlleleFrequencyCalculationModel.Model AFmodel = AlleleFrequencyCalculationModel.Model.EXACT;

    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = DiploidSNPGenotypePriors.HUMAN_HETEROZYGOSITY;

    // control the output
    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confident genotypes (i.e. including ref calls) or just the variants?", required = false)
    public boolean GENOTYPE_MODE = false;

    @Argument(fullName = "output_all_callable_bases", shortName = "all_bases", doc = "Should we output all callable bases?", required = false)
    public boolean ALL_BASES_MODE = false;

    @Argument(fullName = "standard_min_confidence_threshold_for_calling", shortName = "stand_call_conf", doc = "The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be called", required = false)
    public double STANDARD_CONFIDENCE_FOR_CALLING = 30.0;

    @Argument(fullName = "standard_min_confidence_threshold_for_emitting", shortName = "stand_emit_conf", doc = "The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be emitted (and filtered if less than the calling threshold)", required = false)
    public double STANDARD_CONFIDENCE_FOR_EMITTING = 30.0;

    @Argument(fullName = "trigger_min_confidence_threshold_for_calling", shortName = "trig_call_conf", doc = "The minimum phred-scaled confidence threshold at which variants at 'trigger' track sites should be called", required = false)
    public double TRIGGER_CONFIDENCE_FOR_CALLING = 30.0;

    @Argument(fullName = "trigger_min_confidence_threshold_for_emitting", shortName = "trig_emit_conf", doc = "The minimum phred-scaled confidence threshold at which variants at 'trigger' track sites should be emitted (and filtered if less than the calling threshold)", required = false)
    public double TRIGGER_CONFIDENCE_FOR_EMITTING = 30.0;

    @Argument(fullName = "noSLOD", shortName = "nsl", doc = "If provided, we will not calculate the SLOD", required = false)
    public boolean NO_SLOD = false;


    // control the error modes
    @Hidden
    @Argument(fullName = "assume_single_sample_reads", shortName = "single_sample", doc = "The single sample that we should assume is represented in the input bam (and therefore associate with all reads regardless of whether they have read groups)", required = false)
    public String ASSUME_SINGLE_SAMPLE = null;


    // control the various parameters to be used
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public int MIN_BASE_QUALTY_SCORE = 17;

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for calling", required = false)
    public int MIN_MAPPING_QUALTY_SCORE = 20;

    @Argument(fullName = "max_mismatches_in_40bp_window", shortName = "mm40", doc = "Maximum number of mismatches within a 40 bp window (20bp on either side) around the target position for a read to be used for calling", required = false)
    public int MAX_MISMATCHES = 3;

    @Argument(fullName = "use_reads_with_bad_mates", shortName = "bad_mates", doc = "Use reads whose mates are mapped excessively far away for calling", required = false)
    public boolean USE_BADLY_MATED_READS = false;

    @Argument(fullName = "max_deletion_fraction", shortName = "deletions", doc = "Maximum fraction of reads with deletions spanning this locus for it to be callable [to disable, set to < 0 or > 1; default:0.05]", required = false)
    public Double MAX_DELETION_FRACTION = 0.05;

    @Argument(fullName = "get_indel_alleles_from_vcf", shortName = "getIndelAllelesFromVCF", doc = "Get reference/alt alleles for indel genotyping from VCF", required = false)
    public boolean GET_ALLELES_FROM_VCF = false;

    @Argument(fullName = "min_indel_count_for_genotyping", shortName = "minIndelCnt", doc = "Minimum number of consensus indels required to trigger genotyping run", required = false)
    public int MIN_INDEL_COUNT_FOR_GENOTYPING = 5;

    @Argument(fullName = "indel_heterozygosity", shortName = "indelHeterozygosity", doc = "Heterozygosity for indel calling", required = false)
    public double INDEL_HETEROZYGOSITY = 1.0/8000;

    public UnifiedArgumentCollection clone() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();

        uac.GLmodel = GLmodel;
        uac.heterozygosity = heterozygosity;
        uac.GENOTYPE_MODE = GENOTYPE_MODE;
        uac.ALL_BASES_MODE = ALL_BASES_MODE;
        uac.NO_SLOD = NO_SLOD;
        uac.ASSUME_SINGLE_SAMPLE = ASSUME_SINGLE_SAMPLE;
        uac.STANDARD_CONFIDENCE_FOR_CALLING = STANDARD_CONFIDENCE_FOR_CALLING;
        uac.STANDARD_CONFIDENCE_FOR_EMITTING = STANDARD_CONFIDENCE_FOR_EMITTING;
        uac.TRIGGER_CONFIDENCE_FOR_CALLING = TRIGGER_CONFIDENCE_FOR_CALLING;
        uac.TRIGGER_CONFIDENCE_FOR_EMITTING = TRIGGER_CONFIDENCE_FOR_EMITTING;
        uac.MIN_BASE_QUALTY_SCORE = MIN_BASE_QUALTY_SCORE;
        uac.MIN_MAPPING_QUALTY_SCORE = MIN_MAPPING_QUALTY_SCORE;
        uac.MAX_MISMATCHES = MAX_MISMATCHES;
        uac.USE_BADLY_MATED_READS = USE_BADLY_MATED_READS;
        uac.MAX_DELETION_FRACTION = MAX_DELETION_FRACTION;
        uac.GET_ALLELES_FROM_VCF = GET_ALLELES_FROM_VCF;
        uac.MIN_INDEL_COUNT_FOR_GENOTYPING = MIN_INDEL_COUNT_FOR_GENOTYPING;
        uac.INDEL_HETEROZYGOSITY = INDEL_HETEROZYGOSITY;

        return uac;
    }


}
