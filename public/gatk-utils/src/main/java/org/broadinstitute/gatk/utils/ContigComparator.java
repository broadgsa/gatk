/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils;

import htsjdk.samtools.SAMSequenceDictionary;

import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/23/11
 * Time: 6:07 PM
 *
 * Contig comparator -- sorting contigs like Picard
 *
 *   This is very useful if you want to output your text files or manipulate data in the usual chromosome ordering :
 *    1
 *    2
 *    3
 *    ...
 *    21
 *    22
 *    X
 *    Y
 *    GL***
 *    ...
 * Just use this comparator in any SortedSet class constructor and your data will be sorted like in the BAM file.
 */
public class ContigComparator implements Comparator<String> {
    final SAMSequenceDictionary dict;

    public ContigComparator(final SAMSequenceDictionary dict) {
        if ( dict == null ) throw new IllegalArgumentException("dict cannot be null");
        this.dict = dict;
    }

    @Override
    public int compare(final String chr1, final String chr2) {
        final int index1 = getIndex(chr1);
        final int index2 = getIndex(chr2);
        return Integer.valueOf(index1).compareTo(index2);
    }

    /**
     * Convert contig to its index in the dict, or throw an exception if it's not found or is null
     * @param chr the contig
     */
    private int getIndex(final String chr) {
        if ( chr == null ) throw new IllegalArgumentException("chr is null");
        final int index = dict.getSequenceIndex(chr);
        if ( index == -1 ) throw new IllegalArgumentException("Unknown contig " + chr);
        return index;
    }
}
