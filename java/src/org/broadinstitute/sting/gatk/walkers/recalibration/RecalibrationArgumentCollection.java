package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.cmdLine.Argument;

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
 * Date: Nov 27, 2009
 *
 * A collection of the arguments that are common to both CovariateCounterWalker and TableRecalibrationWalker.
 * This set of arguments will also be passed to the constructor of every Covariate when it is instantiated.
 */

public class RecalibrationArgumentCollection {

    //////////////////////////////////
    // Shared Command Line Arguments
    //////////////////////////////////
    @Argument(fullName="recal_file", shortName="recalFile", required=false, doc="Filename for the outputted covariates table recalibration file")
    public String RECAL_FILE = "output.recal_data.csv";
    @Argument(fullName = "use_original_quals", shortName="OQ", doc="If provided, we will use the quals from the original qualities OQ attribute field instead of the quals in the regular QUALS field", required=false)
    public boolean USE_ORIGINAL_QUALS = false;
    @Argument(fullName="default_read_group", shortName="dRG", required=false, doc="If a read has no read group then default to the provided String.")
    public String DEFAULT_READ_GROUP = ReadGroupCovariate.defaultReadGroup;
    @Argument(fullName="default_platform", shortName="dP", required=false, doc="If a read has no platform then default to the provided String. Valid options are illumina, 454, and solid.")
    public String DEFAULT_PLATFORM = "Illumina";
    @Argument(fullName="force_read_group", shortName="fRG", required=false, doc="If provided, the read group ID of EVERY read will be forced to be the provided String. This is useful to collapse all data into a single read group.")
    public String FORCE_READ_GROUP = null;
    @Argument(fullName="force_platform", shortName="fP", required=false, doc="If provided, the platform of EVERY read will be forced to be the provided String. Valid options are illumina, 454, and solid.")
    public String FORCE_PLATFORM = null;
    @Argument(fullName="solid_recal_mode", shortName="sMode", required = false, doc="How should we recalibrate solid bases in which the reference was inserted? EXPERIMENTAL OPTION. USE AT OWN RISK. Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, COUNT_AS_MISMATCH, or REMOVE_REF_BIAS")
    public String SOLID_RECAL_MODE = "DO_NOTHING";
    @Argument(fullName = "window_size_nqs", shortName="nqs", doc="The window size used by MinimumNQSCovariate for its calculation", required=false)
    public int WINDOW_SIZE = 5;
    @Argument(fullName = "homopolymer_nback", shortName="nback", doc="The number of previous bases to look at in HomopolymerCovariate", required=false)
    public int HOMOPOLYMER_NBACK = 8;

    public boolean checkSolidRecalMode() {
        return ( SOLID_RECAL_MODE.equalsIgnoreCase("DO_NOTHING") || SOLID_RECAL_MODE.equalsIgnoreCase("SET_Q_ZERO") ||
                 SOLID_RECAL_MODE.equalsIgnoreCase("SET_Q_ZERO_BASE_N") || SOLID_RECAL_MODE.equalsIgnoreCase("COUNT_AS_MISMATCH") || SOLID_RECAL_MODE.equalsIgnoreCase("REMOVE_REF_BIAS") ); 
    }
}
