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

import static org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel.Model.*;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;

import java.util.Set;


public class GenotypeCalculationModelFactory {
    //private GenotypeCalculationModelFactory() {} // cannot be instantiated

    public static GenotypeCalculationModel.Model getGenotypeCalculationModel(final String name) {
        GenotypeCalculationModel.Model m = valueOf(name);
        if ( m == null )
            throw new RuntimeException("Unexpected GenotypeCalculationModel " + name);
        else
            return m;
    }

    /**
     * General and correct way to create GenotypeCalculationModel objects for arbitrary models
     *
     * @param UAC           The unified argument  collection
     * @param samples       samples in bam
     * @param out           output
     * @return model
     */
    public static GenotypeCalculationModel makeGenotypeCalculation(UnifiedArgumentCollection UAC,
                                                                   Set<String> samples,
                                                                   GenotypeWriter out) {
        switch ( UAC.genotypeModel ) {
            case SINGLE_SAMPLE: return new SingleSampleGenotypeCalculationModel(UAC.baseModel, samples, UAC.defaultPlatform, out, UAC.GENOTYPE, UAC.LOD_THRESHOLD, UAC.MAX_DELETIONS, UAC.VERBOSE);
            //case MULTI_SAMPLE_EM: return new MultiSampleEMGenotypeCalculationModel(UAC.baseModel, samples, UAC.defaultPlatform, out, UAC.GENOTYPE, UAC.LOD_THRESHOLD, UAC.MAX_DELETIONS, UAC.VERBOSE);
            case MULTI_SAMPLE_ALL_MAFS: return new MultiSampleAllMAFsGenotypeCalculationModel(UAC.baseModel, samples, UAC.defaultPlatform, out, UAC.GENOTYPE, UAC.LOD_THRESHOLD, UAC.MAX_DELETIONS, UAC.VERBOSE);
            default: throw new RuntimeException("Unexpected GenotypeCalculationModel " + UAC.genotypeModel);
        }
    }
}