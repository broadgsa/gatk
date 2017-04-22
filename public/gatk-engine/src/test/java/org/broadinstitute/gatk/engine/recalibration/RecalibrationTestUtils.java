/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration;

import org.broadinstitute.gatk.engine.recalibration.covariates.*;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 12/23/12
 * Time: 1:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class RecalibrationTestUtils {
    public static Covariate[] makeInitializedStandardCovariates() {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
        final Covariate[] covariates = new Covariate[4];
        covariates[0] = new ReadGroupCovariate();
        covariates[1] = new QualityScoreCovariate();
        covariates[2] = new ContextCovariate();
        covariates[3] = new CycleCovariate();
        for ( Covariate cov : covariates ) cov.initialize(RAC);
        return covariates;
    }
}
