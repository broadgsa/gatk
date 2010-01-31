/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.cmdLine.Argument;


public class UnifiedArgumentCollection {

    // control the various models to be used
    @Argument(fullName = "genotype_model", shortName = "gm", doc = "Genotype calculation model to employ -- JOINT_ESTIMATE is currently the default, while POOLED and EM_POINT_ESTIMATE are available.", required = false)
    public GenotypeCalculationModel.Model genotypeModel = GenotypeCalculationModel.Model.JOINT_ESTIMATE;

    @Argument(fullName = "base_model", shortName = "bm", doc = "Base substitution model to employ -- EMPIRICAL is the recommended default, but it's possible to select the ONE_STATE and THREE_STATE models for comparison purposes", required = false)
    public BaseMismatchModel baseModel = BaseMismatchModel.EMPIRICAL;

    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = DiploidGenotypePriors.HUMAN_HETEROZYGOSITY;

    @Argument(fullName = "poolSize", shortName = "ps", doc = "Number of individuals in the pool (for POOLED model only)", required = false)
    public int POOLSIZE = 0;

    // control the output
    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confident genotypes (i.e. including ref calls) or just the variants?", required = false)
    public boolean GENOTYPE = false;

    @Argument(fullName = "output_all_callable_bases", shortName = "all_bases", doc = "Should we output all callable bases?", required = false)
    public boolean ALL_BASES = false;

    @Argument(fullName = "noSLOD", shortName = "nsl", doc = "If provided, we will not calculate the SLOD", required = false)
    public boolean NO_SLOD = false;

    @Argument(fullName = "include_experimental_annotations", shortName = "exp", doc = "Annotate calls with all annotations, including experimental ones", required = false)
    public boolean ALL_ANNOTATIONS = false;


    // control the error modes
    @Argument(fullName = "assume_single_sample_reads", shortName = "single_sample", doc = "The single sample that we should assume is represented in the input bam (and therefore associate with all reads regardless of whether they have read groups)", required = false)
    public String ASSUME_SINGLE_SAMPLE = null;

    @Argument(fullName = "platform", shortName = "pl", doc = "Causes the genotyper to assume that reads without PL header TAG are this platform.  Defaults to null, indicating that the system will throw a runtime exception when such reads are detected", required = false)
    public EmpiricalSubstitutionProbabilities.SequencerPlatform defaultPlatform = null;


    // control the various parameters to be used
    @Argument(fullName = "min_confidence_threshold", shortName = "confidence", doc = "The phred-scaled confidence threshold by which variants should be filtered", required = false)
    public double CONFIDENCE_THRESHOLD = 50.0;

    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public int MIN_BASE_QUALTY_SCORE = 10;

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for calling", required = false)
    public int MIN_MAPPING_QUALTY_SCORE = 10;

    @Argument(fullName = "max_mismatches_in_40bp_window", shortName = "mm40", doc = "Maximum number of mismatches within a 40 bp window (20bp on either side) around the target position for a read to be used for calling", required = false)
    public int MAX_MISMATCHES = 3;

    @Argument(fullName = "use_reads_with_bad_mates", shortName = "bad_mates", doc = "Use reads whose mates are mapped excessively far away for calling", required = false)
    public boolean USE_BADLY_MATED_READS = false;

    @Argument(fullName = "max_deletion_fraction", shortName = "deletions", doc = "Maximum fraction of reads with deletions spanning this locus for it to be callable [to disable, set to < 0 or > 1; default:0.05]", required = false)
    public Double MAX_DELETION_FRACTION = 0.05;
}
