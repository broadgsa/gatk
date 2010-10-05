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
import org.apache.log4j.Logger;

import java.util.Set;
import java.io.PrintStream;


public class GenotypeCalculationModelFactory {

    public static GenotypeCalculationModel.Model getGenotypeCalculationModel(final String name) {
        return valueOf(name);
    }

    /**
     * General and correct way to create GenotypeCalculationModel objects for arbitrary models
     *
     * @param samples       samples in bam
     * @param logger        logger
     * @param UAC           the unified argument collection
     * @param verboseWriter verbose writer
     *
     * @return model
     */
    public static GenotypeCalculationModel makeGenotypeCalculation(Set<String> samples,
                                                                   Logger logger,
                                                                   UnifiedArgumentCollection UAC,
                                                                   PrintStream verboseWriter) {
        GenotypeCalculationModel gcm;
        switch ( UAC.genotypeModel ) {
            case JOINT_ESTIMATE:
                gcm = new DiploidGenotypeCalculationModel(samples.size(), UAC.heterozygosity);
                break;
            case DINDEL:
                throw new UnsupportedOperationException("The Dindel-based genotype likelihoods model is not currently supported");
                //gcm = new SimpleIndelCalculationModel();
                //break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + UAC.genotypeModel);
        }

        gcm.initialize(samples, logger, UAC, verboseWriter);
        return gcm;
    }
}