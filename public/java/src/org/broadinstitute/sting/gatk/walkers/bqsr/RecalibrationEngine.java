package org.broadinstitute.sting.gatk.walkers.bqsr;

import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.recalibration.RecalibrationTables;
import org.broadinstitute.sting.utils.recalibration.covariates.Covariate;
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
public interface RecalibrationEngine {
    /**
     * Initialize the recalibration engine
     *
     * Called once before any calls to updateDataForRead are made.  The engine should prepare itself
     * to handle any number of updateDataForRead calls containing ReadRecalibrationInfo containing
     * keys for each of the covariates provided.
     *
     * The engine should collect match and mismatch data into the recalibrationTables data.
     *
     * @param covariates an array of the covariates we'll be using in this engine, order matters
     * @param recalibrationTables the destination recalibrationTables where stats should be collected
     */
    public void initialize(final Covariate[] covariates, final RecalibrationTables recalibrationTables);

    /**
     * Update the recalibration statistics using the information in recalInfo
     * @param recalInfo data structure holding information about the recalibration values for a single read
     */
    @Requires("recalInfo != null")
    public void updateDataForRead(final ReadRecalibrationInfo recalInfo);

    /**
     * Finalize, if appropriate, all derived data in recalibrationTables.
     *
     * Called once after all calls to updateDataForRead have been issued.
     */
    public void finalizeData();
}
