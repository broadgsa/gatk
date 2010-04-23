package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.playground.utils.report.utils.Node;

import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class RFormat
 *
 * a format for outputting R data - experimental
 */
public class RFormat extends TableBasedFormat {
    private final String DIVIDER = ",";
    private static final String extension = ".csv";

    /**
     * format the string according to our internal rules
     *
     * @param str the string to format
     * @return a string, properly formatted
     */
    @Override
    public String formatColumn(String str) {
        return str+DIVIDER;
    }

    /**
     * should we add readability marks?
     *
     * @return true if we should (line breaks, etc)
     */
    @Override
    public boolean addReadabilityMarks() {
        return false;
    }

    /**
     * a string to prepend for header lines
     *
     * @return a string, blank if no string to be appended
     */
    @Override
    public String headerIndicator() {
        return "#";
    }

    /**
     * should we split the seperate files by analysis
     *
     * @return
     */
    @Override
    public boolean splitFilesByAnalysis() {
        return true;
    }

    /**
     * what extension do we want our files to have
     *
     * @return a string of the extension
     */
    @Override
    public String extension() {
        return extension;
    }
}
