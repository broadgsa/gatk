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

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Sep 29, 2009
 * Time: 12:07:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class UnifiedArgumentCollection {

    // control the various models to be used
    @Argument(fullName = "genotypeModel", shortName = "gm", doc = "Genotype calculation model to employ -- SINGLE_SAMPLE and MULTI_SAMPLE_EM are the recommended choices, but it's possible to select MULTI_SAMPLE_ALL_MAFs", required = true)
    public GenotypeCalculationModel.Model genotypeModel = null;

    @Argument(fullName = "baseModel", shortName = "bm", doc = "Base substitution model to employ -- EMPIRICAL is the recommended default, but it's possible to select the ONE_STATE and THREE_STATE models for comparison purposes", required = false)
    public BaseMismatchModel baseModel = BaseMismatchModel.EMPIRICAL;

    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = DiploidGenotypePriors.HUMAN_HETEROZYGOSITY;


    // control the output
    @Argument(fullName = "genotype", shortName = "genotype", doc = "Should we output confident genotypes or just the variants?", required = false)
    public boolean GENOTYPE = false;

    @Argument(fullName = "verbose", shortName = "v", doc = "EXPERIMENTAL", required = false)
    public boolean VERBOSE = false;


    // control the various parameters to be used
    @Argument(fullName = "lod_threshold", shortName = "lod", doc = "The lod threshold on which variants should be filtered", required = false)
    public double LOD_THRESHOLD = Double.MIN_VALUE;

    @Argument(fullName = "platform", shortName = "pl", doc = "Causes the genotyper to assume that reads without PL header TAG are this platform.  Defaults to null, indicating that the system will throw a runtime exception when such reads are detected", required = false)
    public EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform defaultPlatform = null;

    @Argument(fullName = "maxDeletions", shortName = "deletions", doc = "Maximum reads with deletions spanning this locus for it to be callable [default:1]", required = false)
    public Integer MAX_DELETIONS = 1;

    @Argument(fullName = "maxCoverage", shortName = "maxCoverage", doc = "Maximum reads at this locus for it to be callable; to disable, provide value < 1 [default:10,000]", required = false)
    public Integer MAX_READS_IN_PILEUP = 10000;

    //@Argument(fullName = "disableCache", doc = "[ADVANCED] If true, we won't use the caching system.  This argument is for testing purposes only", required = false)
    //public boolean disableCache = false;    

}
