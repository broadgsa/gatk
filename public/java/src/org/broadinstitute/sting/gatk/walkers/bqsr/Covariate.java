package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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
 * The Covariate interface. A Covariate is a feature used in the recalibration that can be picked out of the read.
 * In general most error checking and adjustments to the data are done before the call to the covariates getValue methods in order to speed up the code.
 * This unfortunately muddies the code, but most of these corrections can be done per read while the covariates get called per base, resulting in a big speed up.
 */

public interface Covariate {
    /**
     * Initialize any member variables using the command-line arguments passed to the walker
     *
     * @param RAC the recalibration argument collection
     */
    public void initialize(final RecalibrationArgumentCollection RAC);

    /**
     * Calculates covariate values for all positions in the read.
     *
     * @param read the read to calculate the covariates on.
     * @return all the covariate values for every base in the read.
     */
    public CovariateValues getValues(final GATKSAMRecord read);

    /**
     * Used to get the covariate's value from input csv file during on-the-fly recalibration
     *
     * @param str the key in string type (read from the csv)
     * @return the key in it's correct type.
     */
    public Object getValue(final String str);

    /**
     * Converts the internal representation of the key to String format for file output.
     *
     * @param key the Long representation of the key
     * @return a string representation of the key
     */
    public String formatKey(final Long key);

    /**
     * Converts an Object key into a Long key using only the lowest numberOfBits() bits
     *
     * Only necessary for on-the-fly recalibration when you have the object, but need to store it in memory in long format. For counting covariates
     * the getValues method already returns all values in Long format.
     *
     * @param key the object corresponding to the covariate
     * @return a long representation of the object
     */
    public Long longFromKey(final Object key);

    /**
     * Each covariate should determine how many bits are necessary to encode it's data
     *
     * @return The number of bits used to represent the values of this covariate.
     */
    public int numberOfBits();

}

interface RequiredCovariate extends Covariate {}

interface StandardCovariate extends Covariate {}

interface ExperimentalCovariate extends Covariate {}
