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
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;


public class UnifiedArgumentCollection {

    // control the various models to be used
    @Argument(fullName = "genotype_model", shortName = "gm", doc = "Genotype calculation model to employ -- EM_POINT_ESTIMATE is currently the default, while JOINT_ESTIMATE is under development.  An exception will be thrown if an attempt is made to use EM_POINT_ESTIMATE with pooled data.", required = false)
    public GenotypeCalculationModel.Model genotypeModel = GenotypeCalculationModel.Model.EM_POINT_ESTIMATE;

    @Argument(fullName = "base_model", shortName = "bm", doc = "Base substitution model to employ -- EMPIRICAL is the recommended default, but it's possible to select the ONE_STATE and THREE_STATE models for comparison purposes", required = false)
    public BaseMismatchModel baseModel = BaseMismatchModel.EMPIRICAL;

    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = DiploidGenotypePriors.HUMAN_HETEROZYGOSITY;

    @Argument(fullName = "pooled", shortName = "pooled", doc = "Does the input bam represent pooled data (so that genotypes can't be called)?", required = false)
    public boolean POOLED = false;

    @Argument(fullName = "poolSize", shortName = "ps", doc = "Number of individuals in the pool", required = false)
    public int POOLSIZE = 0;


    // control the output
    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "Format to be used to represent variants; default is VCF", required = false)
    public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;

    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confident genotypes or just the variants?", required = false)
    public boolean GENOTYPE = false;

    @Argument(fullName = "output_all_callable_bases", shortName = "all_bases", doc = "Should we output nonconfident variant or confident ref calls too?", required = false)
    public boolean ALL_BASES = false;

    @Argument(fullName = "verbose_mode", shortName = "verbose", doc = "File to print all of the annotated and detailed debugging output", required = false)
    public String VERBOSE = null;


    // control the error modes
    @Argument(fullName = "assume_single_sample_reads", shortName = "single_sample", doc = "The single sample that we should assume is represented in the input bam (and therefore associate with all reads regardless of whether they have read groups)", required = false)
    public String ASSUME_SINGLE_SAMPLE = null;

    @Argument(fullName = "platform", shortName = "pl", doc = "Causes the genotyper to assume that reads without PL header TAG are this platform.  Defaults to null, indicating that the system will throw a runtime exception when such reads are detected", required = false)
    public EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform defaultPlatform = null;


    // control the various parameters to be used
    @Argument(fullName = "lod_threshold", shortName = "lod", doc = "The lod threshold by which variants should be filtered (for EM_POINT_ESTIMATE model)", required = false)
    public double LOD_THRESHOLD = Double.MIN_VALUE;

    @Argument(fullName = "min_confidence_threshold", shortName = "confidence", doc = "The phred-scaled confidence threshold by which variants should be filtered (for JOINT_ESTIMATE model)", required = false)
    public double CONFIDENCE_THRESHOLD = Double.MIN_VALUE;

    @Argument(fullName = "max_deletions", shortName = "deletions", doc = "Maximum reads with deletions spanning this locus for it to be callable [default:1]", required = false)
    public Integer MAX_DELETIONS = 1;

    @Argument(fullName = "max_coverage", shortName = "coverage", doc = "Maximum reads at this locus for it to be callable; to disable, provide value < 1 [default:10,000]", required = false)
    public Integer MAX_READS_IN_PILEUP = 10000;

    @Argument(fullName = "min_allele_frequency", shortName = "minFreq", doc = "The minimum possible allele frequency in a population (advanced)", required = false)
    public double MINIMUM_ALLELE_FREQUENCY = 1e-8;
}
