/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.illumina;

import java.io.File;
import java.util.Comparator;

/**
 * Comparator for getting Bustard files in "sorted" order for use by the BustardFileParser.  Expected order is
 * by lane in ascending order, then by tile in ascending order, then:
 *      the read 1 qseq file
 *      the read 2 qseq file
 *      the sig2 file
 *
 * IMPORTANT: Currently this class expects to receive ONLY qseq and sig2 files.
 *
 * @author Kathleen Tibbetts
 */
public class BustardFilenameComparator implements Comparator<File> {

    /**
     * Compares its two arguments for order. Returns a negative integer, zero, or a positive integer as
     * the first argument is less than, equal to, or greater than the second.
     *
     * @param file1
     * @param file2
     * @return  a negative integer, zero, or a positive integer as
     *          the first argument is less than, equal to, or greater than the second.
     */
    public int compare(File file1, File file2)
    {
        Integer parts1[] =  parseFileNameParts(file1.getName());
        Integer parts2[] =  parseFileNameParts(file2.getName());

        for (int i = 1; i < parts1.length; i++)
        {
            if (!parts1[i].equals(parts2[i])) {
                return parts1[i].compareTo(parts2[i]);
            }
        }
        return 0;
    }

    /**
     * Utility method that returns an array of integers that represent, in order,
     * the lane, tile, type (0 for qseq files, 1 for sig2 files), and read (if any)
     * represented by the given file name
     *
     * @param name
     * @return  an array of integers that represent, in order,
     *          the lane, tile, type (0 for qseq files, 1 for sig2 files), and read (if any)
     *          represented by the given file name
     */
    private Integer[] parseFileNameParts(String name)
    {
        Integer parts[] = new Integer[4];   // Lane, tile, type, read
        String src[] = name.split("_");
        parts[0] = new Integer(src[1]);     // Lane is always the second part
        if (src[2].length() == 4) {         // Tile is 3rd or fourth
            parts[1] = new Integer(src[2]);
        }
        else {
            parts[1] = new Integer(src[3]);
        }
        parts[2] = (src[src.length-1].equals("qseq.txt")) ? 0 : 1;  // qseq tests are lower
        if (src[2].length() == 1) {  // read is last
            parts[3] = new Integer(src[2]);
        }
        return parts;
    }
}
