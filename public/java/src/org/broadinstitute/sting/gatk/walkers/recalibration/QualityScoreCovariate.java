package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.recalibration.BaseRecalibration;

import java.util.Arrays;

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
 * Date: Nov 3, 2009
 *
 * The Reported Quality Score covariate.
 */

public class QualityScoreCovariate implements RequiredCovariate {

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize( final RecalibrationArgumentCollection RAC ) {
    }

    @Override
    public void getValues( final SAMRecord read, final Comparable[] comparable, final BaseRecalibration.BaseRecalibrationType modelType ) {
        if( modelType == BaseRecalibration.BaseRecalibrationType.BASE_SUBSTITUTION ) {
            byte[] baseQualities = read.getBaseQualities();
            for(int i = 0; i < read.getReadLength(); i++) {
                comparable[i] = (int) baseQualities[i];
            }
        } else { // model == BASE_INSERTION || model == BASE_DELETION
            Arrays.fill(comparable, 45); // Some day in the future when base insertion and base deletion quals exist the samtools API will
                                         // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
        }
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    @Override
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }
}
