package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

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
 * Date: Nov 4, 2009
 *
 * The Mapping Quality covariate.
 */

public class MappingQualityCovariate implements Covariate {

    public MappingQualityCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
    }

    public final Comparable getValue( final ReadHashDatum readDatum, final int offset ) {
    	
    	return readDatum.mappingQuality;
    }
    
    public final Comparable getValue( final String str ) {
        return (int)Integer.parseInt( str ); // cast to primitive int (as opposed to Integer Object) is required so that the return value from the two getValue methods hash to same thing
    }

    public final int estimatedNumberOfBins() {
        return 100;
    }
    
    public String toString() {
        return "Mapping Quality Score";
    }
}
