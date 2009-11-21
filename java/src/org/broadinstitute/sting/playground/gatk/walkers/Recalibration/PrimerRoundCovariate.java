package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.utils.StingException;

import net.sf.samtools.SAMRecord;

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
 * Date: Nov 13, 2009
 *
 * The Primer Round covariate.
 *  For Solexa and 454 this is the same value of the length of the read.
 *  For SOLiD this is different for each position according to http://www3.appliedbiosystems.com/cms/groups/mcb_marketing/documents/generaldocuments/cms_057511.pdf  
 */

public class PrimerRoundCovariate implements Covariate {

	public PrimerRoundCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
    }

    public final Comparable getValue( final ReadHashDatum readDatum, final int offset ) {
        if( readDatum.platform.equalsIgnoreCase( "SOLID" ) ) {
	        return offset % 5; // the primer round according to http://www3.appliedbiosystems.com/cms/groups/mcb_marketing/documents/generaldocuments/cms_057511.pdf
        } else {
            return 1; // nothing to do here because it is always the same
        }

    }

    public final Comparable getValue( final String str ) {
        return (int)Integer.parseInt( str ); // cast to primitive int (as opposed to Integer Object) is required so that the return value from the two getValue methods hash to same thing
    }

    public final int estimatedNumberOfBins() {
        return 5;
    }

    public String toString() {
        return "Primer Round";
    }
}