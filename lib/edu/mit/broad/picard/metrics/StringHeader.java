package edu.mit.broad.picard.metrics;

import edu.mit.broad.picard.util.StringUtil;

/**
 * A simple header who's data type is a single String.  Should not be used for anything other
 * than comments or descriptive text.
 *
 * @author Tim Fennell
 */
public class StringHeader implements Header {
    private String value;

    /** Default constructor. */
    public StringHeader() {}

    /** Constructor that uses the supplied value as the value of the header. */
    public StringHeader(String value) {
        setValue(value);
    }

    public void parse(String in) { value = in.trim(); }
    public String toString() { return value; }

    public String getValue() { return value; }
    public void setValue(String value) { this.value = StringUtil.assertCharactersNotInString(value, '\n'); }

    /** Checks equality on the value of the header. */
    public boolean equals(Object o) {
        if (o != null && o instanceof StringHeader) {
            StringHeader that = (StringHeader) o;
            if (this.value == null) {
                return that.value == null;
            }
            else {
                return this.value.equals(that.value);
            }
        }
        else {
            return false;
        }
    }
}
