package org.broadinstitute.sting.playground.utils.report.templates;

/**
 * the basic comma seperated value format
 */
public class CSVFormat extends TableBasedFormat {
    private final String DIVIDER = ",";

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
     * does the output format want to display line breaks (dotted lines)?
     *
     * @return true if the format uses them
     */
    @Override
    public boolean displayDashedLineBreaks() {
        return false;
    }
}
