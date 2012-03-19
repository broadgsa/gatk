package org.broadinstitute.sting.gatk.walkers.recalibration;

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

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: Mar 22, 2012
 *
 * Object that holds the empirical quality and estimated reported quality values for on-the-fly recalibration. This is a simplification of the RecalDatum object
 */

public class EmpiricalQual {

    private double estimatedQReported; // estimated reported quality score based on combined data's individual q-reporteds and number of observations
    private double empiricalQuality; // the empirical quality for datums that have been collapsed together (by read group and reported quality, for example)

    private EmpiricalQual() {}

    public EmpiricalQual(final double estimatedQReported, final double empiricalQuality) {
        this.estimatedQReported = estimatedQReported;
        this.empiricalQuality = empiricalQuality;
    }

    public final double getEstimatedQReported() {
        return estimatedQReported;
    }

    public final double getEmpiricalQuality() {
        return empiricalQuality;
    }
}