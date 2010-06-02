/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

/**
 * Used to parse files passed to the GenomicAnnotator via the -J arg.
 * The files must be tab-delimited, and the first non-empty/non-commented line
 * must be a header containing column names.
 *
 * More information can be found here: http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator
 */
public class JoinTableParser
{
    public static final String DELIMITER = "\t";

    private List<String> header; //column names parsed out of the header line

    /** Whether to throw an exception if a row contains a different number of columns than the header. */
    private boolean strict;


    /**
     * Constructor.
     *
     * @param source The file to read.
     * @param strict Whether to throw an exception if a row contains a different number of columns than the header.
     */
    public JoinTableParser(boolean strict)  {
        this.strict = strict;
    }

    /**
     * Returns the header and returns it.
     * @param source
     * @return
     */
    public List<String> readHeader(BufferedReader br) throws IOException
    {
        if(header != null) {
            throw new StingException("readHeader(..) called more than once. Header is currently set to: " + header);
        }

        header = Collections.unmodifiableList(parseHeader(br));

        return header;
    }


    /**
     * @return A list containing the column names.
     */
    public List<String> getHeader() {
        return header;
    }


    /**
     * Parses the line into an ArrayList containing the values for each column.
     *
     * @param line
     */
    public ArrayList<String> parseLine(String line) {

        final ArrayList<String> values = Utils.split(line, DELIMITER, header.size());

        if ( strict && values.size() != header.size() ) {
            throw new StingException(String.format("Encountered a row with %d columns which is different from the number or columns in the header: %d\nHeader: " + header + "\nLine: " + values, values.size(), header.size()));
        }

        return values;
    }





    /**
     * Returns the header.
     * @param source The file to read.
     * @return ArrayList containing column names from the header.
     * @throws IOException
     */
    public static ArrayList<String> parseHeader(final BufferedReader br) throws IOException
    {
        ArrayList<String> header = null;

        //find the 1st line that's non-empty and not a comment
        String line = null;
        while( (line = br.readLine()) != null ) {
            line = line.trim();
            if ( line.isEmpty() || line.startsWith("#") ) {
                continue;
            }

            //parse the header
            header = Utils.split(line, DELIMITER);
            break;
        }

        // check that header was found
        if ( header == null ) {
            throw new IllegalArgumentException("No header in " + br + ". All lines are either comments or empty.");
        }

        return header;
    }
}