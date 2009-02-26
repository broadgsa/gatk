package edu.mit.broad.picard.metrics;

/**
 * A header for a metrics file.  A header simply consists of a type and some arbitrary
 * data, but must be able to turn itself into a String and parse it's data back out
 * of that String at a later date.
 *
 * @author Tim Fennell
 */
public interface Header {
    /** Converts the header to a String for persisting to a file. */
    public String toString();

    /** Parses the data contained in the String version of the header. */
    public void parse(String in);

}
