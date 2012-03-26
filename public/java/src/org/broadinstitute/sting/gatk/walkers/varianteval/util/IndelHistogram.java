/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;

/**
 * Simple utility for histogramming indel lengths
 *
 * Based on code from chartl
 *
 * @author Mark DePristo
 * @since 3/21/12
 */
public class IndelHistogram extends TableType {
    private final boolean asFrequencies;
    int nIndels = 0, nTooBigDeletions = 0, nTooBigInsertions = 0;
    private final Integer[] rowKeys;

    private Map<Integer, Double> frequencies = null;
    private final Map<Integer, Integer> counts = new HashMap<Integer, Integer>();
    
    public IndelHistogram(int maxSize, boolean asFrequencies) {
        this.asFrequencies = asFrequencies;
        initializeCounts(maxSize);
        this.rowKeys = new ArrayList<Integer>(counts.keySet()).toArray(new Integer[maxSize]);
    }

    private void initializeCounts(int size) {
        for ( int i = -size; i <= size; i++ ) {
            if ( i != 0 ) counts.put(i, 0);
        }
    }

    @Override
    public String getRowName() {
        return "Length";
    }

    @Override
    public Object[] getColumnKeys() {
        return new String[]{"Count"};
    }

    @Override
    public Object[] getRowKeys() {
        return rowKeys;
    }

    @Override
    public Object getCell(int row, int col) {
        final int key = (Integer)getRowKeys()[row];
        if ( asFrequencies ) {
            if ( frequencies == null ) {
                frequencies = new HashMap<Integer, Double>();
                for ( final int len : counts.keySet() ) {
                    final double value = nIndels == 0 ? 0.0 : counts.get(len) / (1.0 * nIndels);
                    frequencies.put(len, value);
                }
            }
            return frequencies.get(key);
        }
        return counts.get(key);
    }

    public int getnTooBigDeletions() {
        return nTooBigDeletions;
    }

    public int getnTooBigInsertions() {
        return nTooBigInsertions;
    }

    public void update(final Allele ref, final Allele alt) {
        final int alleleSize = alt.length() - ref.length();
        update(alleleSize);
    }
    
    public void update(int len) {
        if ( counts.containsKey(len) ) {
            nIndels++;
            counts.put(len, counts.get(len) + 1);
        } else if ( len < 0 ) {
            nTooBigDeletions++;
        } else {
            nTooBigInsertions++;
        }
    }
}
