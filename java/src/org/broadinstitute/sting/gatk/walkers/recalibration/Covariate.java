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
 * User: rpoplin
 * Date: Oct 30, 2009
 *
 * The Covariate interface. A Covariate is a feature used in the recalibration that can be picked out of the read, offset, and corresponding reference bases
 * In general most error checking and adjustments to the data are done before the call to the covariates getValue methods in order to speed up the code.
 * This unfortunately muddies the code, but most of these corrections can be done per read while the covariates get called per base, resulting in a big speed up.
 */

public interface Covariate {
    public void initialize( RecalibrationArgumentCollection RAC ); // Initialize any member variables using the command-line arguments passed to the walkers
    public Comparable getValue( ReadHashDatum readDatum, int offset ); // Used to pick out the covariate's value from attributes of the read
    public Comparable getValue( String str ); // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public int estimatedNumberOfBins(); // Used to estimate the amount space required for the full data HashMap
}
