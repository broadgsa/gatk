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
