/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 27, 2009
 *
 * A collection of the arguments that are common to both CovariateCounterWalker and TableRecalibrationWalker.
 * This set of arguments will also be passed to the constructor of every Covariate when it is instantiated.
 */

public class RecalibrationArgumentCollection {

    //////////////////////////////////
    // Shared Command Line Arguments
    //////////////////////////////////
    @Hidden
    @Argument(fullName="default_read_group", shortName="dRG", required=false, doc="If a read has no read group then default to the provided String.")
    public String DEFAULT_READ_GROUP = null;
    @Hidden
    @Argument(fullName="default_platform", shortName="dP", required=false, doc="If a read has no platform then default to the provided String. Valid options are illumina, 454, and solid.")
    public String DEFAULT_PLATFORM = null;
    @Hidden
    @Argument(fullName="force_read_group", shortName="fRG", required=false, doc="If provided, the read group ID of EVERY read will be forced to be the provided String. This is useful to collapse all data into a single read group.")
    public String FORCE_READ_GROUP = null;
    @Hidden
    @Argument(fullName="force_platform", shortName="fP", required=false, doc="If provided, the platform of EVERY read will be forced to be the provided String. Valid options are illumina, 454, and solid.")
    public String FORCE_PLATFORM = null;
    @Hidden
    @Argument(fullName = "window_size_nqs", shortName="nqs", doc="The window size used by MinimumNQSCovariate for its calculation", required=false)
    public int WINDOW_SIZE = 5;

    /**
     * This window size tells the module in how big of a neighborhood around the current base it should look for the minimum base quality score.
     */
    @Hidden
    @Argument(fullName = "homopolymer_nback", shortName="nback", doc="The number of previous bases to look at in HomopolymerCovariate", required=false)
    public int HOMOPOLYMER_NBACK = 7;
    @Hidden
    @Argument(fullName = "exception_if_no_tile", shortName="throwTileException", doc="If provided, TileCovariate will throw an exception when no tile can be found. The default behavior is to use tile = -1", required=false)
    public boolean EXCEPTION_IF_NO_TILE = false;

    /**
     * CountCovariates and TableRecalibration accept a --solid_recal_mode <MODE> flag which governs how the recalibrator handles the
     * reads which have had the reference inserted because of color space inconsistencies.
     */
    @Argument(fullName="solid_recal_mode", shortName="sMode", required = false, doc="How should we recalibrate solid bases in which the reference was inserted? Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS")
    public RecalDataManager.SOLID_RECAL_MODE SOLID_RECAL_MODE = RecalDataManager.SOLID_RECAL_MODE.SET_Q_ZERO;

    /**
     * CountCovariates and TableRecalibration accept a --solid_nocall_strategy <MODE> flag which governs how the recalibrator handles
     * no calls in the color space tag. Unfortunately because of the reference inserted bases mentioned above, reads with no calls in
     * their color space tag can not be recalibrated.
     */
    @Argument(fullName = "solid_nocall_strategy", shortName="solid_nocall_strategy", doc="Defines the behavior of the recalibrator when it encounters no calls in the color space. Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ", required=false)
    public RecalDataManager.SOLID_NOCALL_STRATEGY SOLID_NOCALL_STRATEGY = RecalDataManager.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION;
}
