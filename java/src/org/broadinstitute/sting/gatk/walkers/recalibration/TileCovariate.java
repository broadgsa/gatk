/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.StingException;
import net.sf.samtools.SAMRecord;
import net.sf.picard.util.IlluminaUtil;

/**
 * @author alecw@broadinstitute.org
 */

public class TileCovariate implements ExperimentalCovariate {

    private boolean exceptionWhenNoTile = false;

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        exceptionWhenNoTile = RAC.EXCEPTION_IF_NO_TILE;
    }


    // Used to pick out the covariate's value from attributes of the read
    public Comparable getValue(final SAMRecord read, final int offset) {
        Integer tile = IlluminaUtil.getTileFromReadName(read.getReadName());
        if (tile == null) {
            if( exceptionWhenNoTile ) {
                throw new StingException( "Tile number not defined for read: " + read.getReadName() );
            } else {
                return -1;
            }
        }
        return tile;
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public Comparable getValue(final String str) {
        return Integer.parseInt( str );
    }

    public void getValues(SAMRecord read, Comparable[] comparable) {
        final Comparable tileCovariateValue = getValue(read, 0);
        for(int i = 0; i < read.getReadLength(); i++) {
            comparable[i] = tileCovariateValue;
        }
    }
}
