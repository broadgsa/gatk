package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;
import java.util.HashMap;

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
 * The Read Group covariate.
 */

public class ReadGroupCovariate implements RequiredCovariate {
    
    private final HashMap<String, Short> readGroupLookupTable = new HashMap<String, Short>();
    private final HashMap<Short, String> readGroupReverseLookupTable = new HashMap<Short, String>();
    private short nextId = 0;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
    }

    @Override
    public CovariateValues getValues(final GATKSAMRecord read) {
        final int l = read.getReadLength();
        final String readGroupId = read.getReadGroup().getReadGroupId();
        short shortId;
        if (readGroupLookupTable.containsKey(readGroupId)) 
            shortId = readGroupLookupTable.get(readGroupId);
        else {
            shortId = nextId;
            readGroupLookupTable.put(readGroupId, nextId);
            readGroupReverseLookupTable.put(nextId, readGroupId);
            nextId++;
        }
        Short [] readGroups = new Short[l];
        Arrays.fill(readGroups, shortId);
        return new CovariateValues(readGroups, readGroups, readGroups);
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return str;
    }
    
    public final String decodeReadGroup(final short id) {
        return readGroupReverseLookupTable.get(id);
    }
}


