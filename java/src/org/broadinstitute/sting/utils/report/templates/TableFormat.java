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

package org.broadinstitute.sting.utils.report.templates;


/**
 * 
 * @author aaron 
 * 
 * Class TableFormat
 *
 * implements the table (human readable) format.
 */
public class TableFormat extends TableBasedFormat {
    private static final int COLUMN_WIDTH = 25;
    private static final String TBL = "tbl";

    /**
     * format the string according to our internal rules
     *
     * @param str the string to format
     * @return a string, properly formatted
     */
    @Override
    public String formatColumn(String str)  {
        return String.format("%-"+COLUMN_WIDTH+"s",str);
    }

    /**
     * should we add readability marks?
     *
     * @return true if we should (line breaks, etc)
     */
    @Override
    public boolean addReadabilityMarks() {
        return true;
    }

    /**
     * a string to prepend for header lines
     *
     * @return a string, blank if no string to be appended
     */
    @Override
    public String headerIndicator() {
        return "";
    }

    /**
     * should we split the seperate files by analysis
     *
     * @return
     */
    @Override
    public boolean splitFilesByAnalysis() {
        return false;
    }

    /**
     * what extension do we want our files to have
     *
     * @return a string of the extension
     */
    @Override
    public String extension() {
        return TBL;
    }
}

