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

package org.broadinstitute.gatk.tools.walkers.annotator;

import org.broadinstitute.gatk.utils.exceptions.GATKException;

import java.util.Arrays;

/**
 * Created by gauthier on 11/1/16.
 */
public class Histogram {
    private Double binSize;
    private String precisionFormat;
    private String printDelim;
    final private Double BIN_EPSILON = 0.01;

    private CompressedDataList<Integer> dataList = new CompressedDataList<>();

    public Histogram() {
        this.binSize = 0.1;
        precisionFormat = "%.1f";
    }

    public Histogram(final Double binSize) {
        this.binSize = binSize;
        precisionFormat = "%." + Math.round(-Math.log10(binSize)) + "f";
    }

    public void add(final Double d) {
        if (d.isNaN())
            return;
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey))
            dataList.add((int) binKey);
        else
            throw new GATKException("Histogram values are suspiciously extreme.  Failed to add " + d + " to the Histogram.");
    }

   public void add(final Double d, final int count) {
       if (count < 1)
           throw new GATKException("Cannot add non-positive counts to Histogram.");
       long binKey = getBinnedValue(d);
       if (isValidBinKey(binKey))
           dataList.add((int)binKey, count);
       else
           throw new GATKException("Histogram values are suspiciously extreme.  Failed to add " + d + " to the Histogram.");
    }

    public void add(final Histogram h) {
        if (!this.binSize.equals(h.binSize))
            throw new GATKException("Histogram bin sizes are mismatched -- cannot add bin size " + this.binSize + " to " + h.binSize);
        this.dataList.add(h.dataList);
    }

    public Integer get(final Double d) {
        long binKey = getBinnedValue(d);
        if (isValidBinKey(binKey))
            return dataList.getValueCounts().get((int)binKey);
        else
            throw new GATKException("Requested value is suspiciously extreme.  Failed to retrieve " + d + " from the Histogram.");
    }

    /**
     *
     * @return may be null if Histogram is empty
     */

    public Double median() {
        int numItems = 0;
        for(final int count : dataList.valueCounts.values()) {
            numItems += count;
        }
        boolean oddNumberValues = true;
        if(numItems % 2 == 0)
            oddNumberValues = false;
        int medianIndex = (numItems+1)/2;

        int counter = 0;
        Double firstMedian = null;
        for(final Integer key : dataList.valueCounts.keySet()) {
            counter += dataList.valueCounts.get(key);
            if( counter > medianIndex) {
                if (firstMedian == null)
                    return key*binSize;
                else {
                    return (firstMedian+key)/2.0*binSize;
                }
            }
            if( counter == medianIndex) {
                if (oddNumberValues)
                    return key*binSize;
                else {
                    firstMedian = (double) key;
                }
            }
        }
        return null;
    }

    private long getBinnedValue(double d) {
        return Math.round(Math.floor((d+BIN_EPSILON*binSize)/binSize)); //add a little epsilon before division so values exactly on bin boundaries will stay in the same bin
    }

    private boolean isValidBinKey(long binnedValue) {
        return binnedValue <= Integer.MAX_VALUE && binnedValue >= Integer.MIN_VALUE;
    }

    @Override
    public String toString(){
        printDelim = ",";
        String str = "";
        Object[] keys = dataList.valueCounts.keySet().toArray();
        if (keys.length == 0)
            return Double.toString(Double.NaN);
        Arrays.sort(keys);
        for (Object i: keys){
            if(!str.isEmpty())
                str+=printDelim;
            str+=(String.format(precisionFormat,(double)(int)i*binSize)+printDelim+dataList.valueCounts.get(i));  //key i needs to be output with specified precision
        }
        return str;
    }

    public boolean isEmpty() {
        return dataList.isEmpty();
    }
}
