package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.BitSetUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;
import java.util.BitSet;
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
        final String readGroupId = readGroupValueFromRG(read.getReadGroup());
        BitSet rg = bitSetForReadGroup(readGroupId);                            // All objects must output a BitSet, so we convert the "compressed" representation of the Read Group into a bitset
        BitSet[] readGroups = new BitSet[l];
        Arrays.fill(readGroups, rg);
        return new CovariateValues(readGroups, readGroups, readGroups);
    }

    @Override
    public final Object getValue(final String str) {
        return str;
    }

    @Override
    public String keyFromBitSet(BitSet key) {
        return decodeReadGroup((short) BitSetUtils.longFrom(key));
    }

    @Override
    public BitSet bitSetFromKey(Object key) {
        return bitSetForReadGroup((String) key);
    }

    @Override
    public int numberOfBits() {
        return BitSetUtils.numberOfBitsToRepresent(Short.MAX_VALUE);
    }

    private String decodeReadGroup(final short id) {
        return readGroupReverseLookupTable.get(id);
    }

    private BitSet bitSetForReadGroup(String readGroupId) {
        short shortId;
        if (readGroupLookupTable.containsKey(readGroupId))
            shortId = readGroupLookupTable.get(readGroupId);
        else {
            shortId = nextId;
            readGroupLookupTable.put(readGroupId, nextId);
            readGroupReverseLookupTable.put(nextId, readGroupId);
            nextId++;
        }        
        return BitSetUtils.bitSetFrom(shortId);
    }

    /**
     * If the sample has a PU tag annotation, return that. If not, return the read group id.
     *
     * @param rg the read group record
     * @return platform unit or readgroup id
     */
    private String readGroupValueFromRG(GATKSAMReadGroupRecord rg) {
        String platformUnit = rg.getPlatformUnit();
        return platformUnit == null ? rg.getId() : platformUnit;
    }
    
}


