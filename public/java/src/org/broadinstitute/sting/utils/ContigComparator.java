/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils;

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
    private Set<String> specialChrs;

    public ContigComparator() {
        specialChrs = new TreeSet<String>();
        specialChrs.add("X");
        specialChrs.add("Y");
    }

    public int compare(String chr1, String chr2) {
        if (chr1.equals(chr2))
            return 0;

        Integer x = convertStringWithoutException(chr1);
        Integer y = convertStringWithoutException(chr2);
        // both contigs are numbered
        if (x != null && y != null)
            return (x < y) ? -1:1;

        // both contigs are named
        if (x == null && y == null) {
            // both contigs are special contigs or neither contig is a special contigs
            if (specialChrs.contains(chr1) && specialChrs.contains(chr2) || (!specialChrs.contains(chr1) && !specialChrs.contains(chr2)))
                return chr1.compareTo(chr2);
            // one contig is a special and the other is not special
            if (specialChrs.contains(chr1))
                return -1;
            return 1;
        }

        // one contig is named the other is numbered
        if (x != null)
            return -1;
        return 1;
    }

    private Integer convertStringWithoutException(String contig) {
        Integer x = null;
        try {
            x = Integer.decode(contig);
        } catch (NumberFormatException n){}
        return x;
    }
}
