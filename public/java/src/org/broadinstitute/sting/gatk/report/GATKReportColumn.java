/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.report;

import org.apache.commons.lang.math.NumberUtils;

import java.util.*;

/**
 * Holds values for a column in a GATK report table
 */
public class GATKReportColumn extends LinkedHashMap<Object, Object> {
    final private String columnName;
    final private Object defaultValue;
    final private String format;
    final private boolean display;
    final private GATKReportDataType dataType;

    /**
     * Construct the column object, specifying the column name, default value, whether or not the column should be
     * displayed, and the format string. This cannot be null.
     *
     * @param columnName   the name of the column
     * @param defaultValue the default value of the column
     * @param display      if true, the column will be displayed in the final output
     * @param format       format string
     */
    public GATKReportColumn(String columnName, Object defaultValue, boolean display, String format) {
        this.columnName = columnName;
        this.display = display;
        if ( format.equals("") ) {
            this.format = "%s";
            this.dataType = GATKReportDataType.Unknown;
            if ( defaultValue != null ) {
                this.defaultValue = defaultValue;
                //this.dataType = GATKReportDataType.fromObject(defaultValue);
            }
            else {
                this.defaultValue = "";
                //this.dataType = GATKReportDataType.Unknown;
            }
        }
        else {
            this.format = format;
            this.dataType = GATKReportDataType.fromFormatString(format);
            if ( defaultValue == null ) {
                this.defaultValue = dataType.getDefaultValue();
            }
            else {
                this.defaultValue = defaultValue;
            }
        }
    }

    /**
     * Initialize an element in the column with a default value
     *
     * @param primaryKey the primary key position in the column that should be set
     */
    public void initialize(Object primaryKey) {
        this.put(primaryKey, defaultValue);
    }

    /**
     * Return an object from the column, but if it doesn't exist, return the default value.  This is useful when writing
     * tables, as the table gets written properly without having to waste storage for the unset elements (usually the zero
     * values) in the table.
     *
     * @param primaryKey the primary key position in the column that should be retrieved
     * @return the value at the specified position in the column, or the default value if the element is not set
     */
    private Object getWithoutSideEffects(Object primaryKey) {
        if (!this.containsKey(primaryKey)) {
            return defaultValue;
        }

        return this.get(primaryKey);
    }

    /**
     * Return an object from the column, but if it doesn't exist, return the default value.
     *
     * @param primaryKey the primary key position in the column that should be retrieved
     * @return the string value at the specified position in the column, or the default value if the element is not set
     */
    public String getStringValue(Object primaryKey) {
        return formatValue(getWithoutSideEffects(primaryKey));
    }

    /**
     * Return the displayable property of the column.  If true, the column will be displayed in the final output.
     * If not, printing will be suppressed for the contents of the table.
     *
     * @return true if the column will be displayed, false if otherwise
     */
    public boolean isDisplayable() {
        return display;
    }

    /**
     * Get the display width for this column.  This allows the entire column to be displayed with the appropriate, fixed width.
     *
     * @return the format string for this column
     */
    public GATKReportColumnFormat getColumnFormat() {
        int maxWidth = columnName.length();
        GATKReportColumnFormat.Alignment alignment = GATKReportColumnFormat.Alignment.RIGHT;

        for (Object obj : this.values()) {
            if (obj != null) {
                String formatted = formatValue(obj);

                int width = formatted.length();
                if (width > maxWidth) {
                    maxWidth = width;
                }

                if (alignment == GATKReportColumnFormat.Alignment.RIGHT) {
                    if (!isRightAlign(formatted)) {
                        alignment = GATKReportColumnFormat.Alignment.LEFT;
                    }
                }
            }
        }

        return new GATKReportColumnFormat(maxWidth, alignment);
    }

    private static final Collection<String> RIGHT_ALIGN_STRINGS = Arrays.asList(
            "null",
            "NA",
            String.valueOf(Double.POSITIVE_INFINITY),
            String.valueOf(Double.NEGATIVE_INFINITY),
            String.valueOf(Double.NaN));

    /**
     * Check if the value can be right aligned. Does not trim the values before checking if numeric since it assumes
     * the spaces mean that the value is already padded.
     *
     * @param value to check
     * @return true if the value is a right alignable
     */
    protected static boolean isRightAlign(String value) {
        return value == null || RIGHT_ALIGN_STRINGS.contains(value) || NumberUtils.isNumber(value);
    }

    /**
     * Returns a string version of the values.
     *
     * @param obj The object to convert to a string
     * @return The string representation of the column
     */
    private String formatValue(Object obj) {
        String value;
        if (obj == null) {
            value = "null";
        } else if ( dataType.equals(GATKReportDataType.Unknown) &&
                (obj instanceof Double || obj instanceof Float) ) {
            value = String.format("%.8f", obj);
        } else
            value = String.format(format, obj);

        return value;
    }

    public GATKReportDataType getDataType() {
        return dataType;
    }

    public boolean isSameFormat(GATKReportColumn that) {
        return (dataType.equals(that.dataType) &&
                columnName.equals(that.columnName) &&
                display == that.display &&
                format.equals(that.format) &&
                defaultValue.equals(that.defaultValue) );
    }

    boolean equals(GATKReportColumn that) {
        if ( !this.keySet().equals(that.keySet()) ) {
            return false;
        }

        for (Object key : keySet()) {
            Object ValueA = this.get(key);
            Object ValueB = that.get(key);

            //if the value is not equal, (use data type to get the right comparison)
            if (!dataType.isEqual(ValueA, ValueB)) {
                return false;
            }
        }

        return true;
    }

    public String getColumnName() {
        return columnName;
    }

    public String getFormat() {
        if ( dataType.equals(GATKReportDataType.Unknown) ) {
            return "";
        }
        else
            return format;
    }
}
