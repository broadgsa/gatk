/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.text.XReadLines;

/**
 * <pre>
 * This slightly modified TabularROD format is used as input to the
 * GenomicAnnotator.
 * The main differences from TabularROD are:
 * - the delimiter is \t instead of \s+
 * - incomplete records are allowed (eg. where some values are left blank)
 * - the header column is just the first non-comment, non-empty row in the file.
 *    It no longer has to start with the "HEADER" keyword.
 * </pre>
 *
 * More details can be found here:  http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator
 */
public class AnnotatorROD extends TabularROD {

    private static Logger logger = Logger.getLogger(AnnotatorROD.class);

    /** Special column names */
    //public static final String CHRPOS_COLUMN = "chrpos";
    public static final String HAPLOTYPE_REFERENCE_COLUMN = "haplotypeReference";
    public static final String HAPLOTYPE_ALTERNATE_COLUMN = "haplotypeAlternate";
    public static final String HAPLOTYPE_STRAND_COLUMN = "haplotypeStrand";


    private static int parsedRecords = 0;

   /**
    * Constructor.
    *
    * @param name The binding name provided as the first in the list of -B args.
    */
    public AnnotatorROD(String name) {
        super(name, new ArrayList<String>());

        setDelimiter("\t", "\\t");
    }


    /**
     * Walks through the source files looking for the header line, which it
     * returns as a list of strings.
     *
     * @param source
     * @return
     */
    public Object initialize(final File source) throws FileNotFoundException {
        ArrayList<String> header = AnnotatorROD.readHeader(source);
        return header;
    }




    /**
     * Finds and parses the header string in the file.
     * @param source The source file.
     * @return A List of the column names parsed out of the header.
     * @throws FileNotFoundException
     */
    public static ArrayList<String> readHeader(final File source) throws FileNotFoundException {
        ArrayList<String> header = null;
        XReadLines reader = new XReadLines(source);

        //find the 1st line that's non-empty and not a comment
        for ( String line : reader ) {
            line = line.trim();
            if ( line.isEmpty() || line.startsWith("#") ) {
                continue;
            }

            header = new ArrayList<String>(Arrays.asList(line.split(DELIMITER_REGEX)));
            break;
        }

        // check that we found the header
        if ( header == null ) {
            throw new RuntimeException("No header in " + source + ". All lines are either comments or empty.");
        }

        logger.info(String.format("Found header line containing %d columns:\n[%s]", header.size(), Utils.join("\t", header)));

        try {
            reader.close();
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }

        return header;
    }


    /**
     * Used by ROD management system to set the data in this ROD associated with a line in a rod
     *
     * @param headerObj
     * @param parts
     * @return
     * @throws IOException
     */
    public boolean parseLine(final Object headerObj, final String[] parts) throws IOException {
        ArrayList<String> header = (ArrayList<String>)(headerObj);

        //save the header in the super-class
        getHeader().addAll(header);

        if ( parts.length == 0 || parts[0].startsWith("#") || header.get(0).equals(parts[0]) /* Skip the header line */ )
            return false;


        if ( header.size() < parts.length) {
            throw new IOException(String.format("Encountered line with more columns than have names in the header. Header has %d columns, this line has %d columns.", header.size(), parts.length));
        }

        for ( int i = 0; i < parts.length; i++ ) {
            put(header.get(i), parts[i]);
        }

        if ( printRecordsParsed ) System.out.printf("Parsed %d records %s%n", ++parsedRecords, this);

        return true;
    }
}
