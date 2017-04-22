/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration.covariates;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;


public class RepeatUnitAndLengthCovariate extends RepeatCovariate {

    @Requires({"repeatLength>=0", "repeatFromUnitAndLength != null"})
    @Ensures("result != null")
    protected String getCovariateValueFromUnitAndLength(final byte[] repeatFromUnitAndLength, final int repeatLength) {
        return new String(repeatFromUnitAndLength) + String.format("%d",repeatLength);
    }

    @Override
    public synchronized int maximumKeyValue() {
        // Synchronized so that we don't query table size while the tables are being updated
        //return repeatLookupTable.size() - 1;
        // max possible values of covariate: for repeat unit, length is up to MAX_STR_UNIT_LENGTH,
        // so we have 4^MAX_STR_UNIT_LENGTH * MAX_REPEAT_LENGTH possible values
        return (1<<(2*MAX_STR_UNIT_LENGTH)) * MAX_REPEAT_LENGTH +1;
    }

}
