package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.playground.utils.report.utils.Node;
import org.broadinstitute.sting.utils.Pair;

import java.io.*;
import java.util.*;


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
     * does the output format want to display line breaks (dotted lines)?
     *
     * @return true if the format uses them
     */
    @Override
    public boolean displayDashedLineBreaks() {
        return true;
    }
}

