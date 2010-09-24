package org.broadinstitute.sting.utils.report.templates;

/**
 * the basic comma separated value format
 */
public class CSVFormat extends TableBasedFormat {
    private static final String DIVIDER = ",";
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
     * should we split the separate files by analysis
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
        return extension;
    }
}
