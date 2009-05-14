package org.broadinstitute.sting.gatk.refdata;

import java.util.*;
import java.util.regex.MatchResult;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.xReadLines;
import org.apache.log4j.Logger;

/**
 * Class for representing arbitrary reference ordered data sets
 *
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 *
 * System for interacting with tabular formatted data of the following format:
 *
 * # comment line
 * # must include HEADER KEYWORD
 * HEADER COL1 ... COLN
 * chr:pos data1 ... dataN
 *
 * The system supports the rod interface.  You can just access tabularRODs through the normal ROD system.
 *
 * You can also write your own files, as such:
 *
 * ArrayList<String> header = new ArrayList<String>(Arrays.asList("HEADER", "col1", "col2", "col3"));
 * assertTrue(TabularROD.headerString(header).equals("HEADER\tcol1\tcol2\tcol3"));
 * String rowData = String.format("%d %d %d", 1, 2, 3);
 * TabularROD row = new TabularROD("myName", header, new GenomeLoc("chrM", 1), rowData.split(" "));
 * assertTrue(row.toString().equals("chrM:1\t1\t2\t3"));
 */
public class TabularROD extends BasicReferenceOrderedDatum implements Map<String, String> {
    private static Logger logger = Logger.getLogger(TabularROD.class);

    private GenomeLoc loc;
    private HashMap<String, String> attributes;
    private ArrayList<String> header;

    public static String DEFAULT_DELIMITER = "\t";
    public static String DEFAULT_DELIMITER_REGEX = "\\s+";

    public static String DELIMITER = DEFAULT_DELIMITER;
    public static String DELIMITER_REGEX = DEFAULT_DELIMITER_REGEX;

    private static int MAX_LINES_TO_LOOK_FOR_HEADER = 1000;
    private static Pattern HEADER_PATTERN = Pattern.compile("^\\s*HEADER.*");
    private static Pattern COMMENT_PATTERN = Pattern.compile("^#.*");

    /**
     * Set the global tabular ROD delimiter and the regex to split the delimiter.
     *
     * The delimiter to put between fields, while the regex is used to split lines
     * 
     * @param delimiter
     * @param delimeterRegex
     */
    public static void setDelimiter(final String delimiter, final String delimeterRegex) {
        DELIMITER = delimiter;
        DELIMITER_REGEX = delimeterRegex;
    }

    /**
     * Returns a parsable string representation for the
     * @param header
     */
    public static String headerString(final ArrayList<String> header) {
        requireGoodHeader(header);
        return Utils.join(DELIMITER, header);
    }

    /**
     * Returns a comment line containing the *single line* string msg
     * 
     * @param msg
     * @return
     */
    public static String commentString(final String msg) {
        return "# " + msg;
    }

    private static boolean headerIsGood(final ArrayList<String> header) {
        if ( header.size() == 0 ) return false;
        if ( ! header.get(0).equals("HEADER") ) return false;
        return true;
    }

    private static void requireGoodHeader(final ArrayList<String> header) {
        if ( ! headerIsGood(header) )
            throw new RuntimeException("Header must begin with HEADER keyword");
    }

    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public TabularROD(final String name) {
        super(name);
        attributes = new HashMap<String, String>();
    }

    /**
     * Make a new TabularROD with name, using header columns header, at loc, without any bound data.  Data
     * must be bound to each corresponding header[i] field before the object is really usable.
     *
     * @param name
     * @param header
     * @param loc
     */
    public TabularROD(final String name, ArrayList<String> header, GenomeLoc loc) {
        this(name);
        this.header = header;
        this.loc = loc;
        requireGoodHeader(this.header);
    }

    /**
     * Make a new TabularROD with name, using header columns header, at loc, with data associated with the
     * header columns.  data and header are assumed to be in the same order, and bindings will be established
     * from header[i] = data[i].  The TabularROD at this stage can be printed, manipulated, it is considered
     * a full fledged, initialized object.
     *
     * @param name
     * @param header
     * @param loc
     * @param data
     */
    public TabularROD(final String name, ArrayList<String> header, GenomeLoc loc, String[] data) {
        this(name, header, loc);
        
        if ( header.size() != data.length + 1 )
            throw new RuntimeException(String.format("Incorrect tabular data format: header has %d columns but %d data elements were provided: %s",
                                                    header.size(), data.length, Utils.join("\t", data)));
        for ( int i = 0; i < data.length; i++ ) {
            put(header.get(i+1), data[i]);
        }
    }

    /**
     * Walks through the source files looking for the header line, which it returns as a
     * list of strings.
     * 
     * @param source
     * @return
     */
    public Object initialize(final File source) throws FileNotFoundException {
        List<String> header = null;
        int linesLookedAt = 0;
        xReadLines reader = new xReadLines(source);

        for ( String line : reader ) {
            Matcher m = HEADER_PATTERN.matcher(line);
            if ( m.matches() ) {
                //System.out.printf("Found a header line: %s%n", line);
                header = new ArrayList<String>(Arrays.asList(line.split(DELIMITER_REGEX)));
                //System.out.printf("HEADER IS %s%n", Utils.join(":", header));
            }

            if ( linesLookedAt++ > MAX_LINES_TO_LOOK_FOR_HEADER )
                break;
        }

        // check that we found the header
        if ( header != null ) {
            logger.debug(String.format("HEADER IS %s%n", Utils.join(":", header)));
        } else {
            throw new RuntimeException("Couldn't find header line in TabularROD!");
        }

        return header;
    }

    // ----------------------------------------------------------------------
    //
    // ROD accessors
    //
    // ----------------------------------------------------------------------
    public GenomeLoc getLocation() {
        return loc;
    }

    public ArrayList<String> getHeader() {
        return header;
    }

    public String get(final Object key) {
        return attributes.get(key);
    }

    public String put(final String key, final String object) {
        return attributes.put(key, object);
    }

    public boolean containsKey(final Object key) {
        return attributes.containsKey(key);
    }

    public HashMap<String,String> getAttributes() {
        return attributes;
    }

    public String getAttributeString() {
        List<String> strings = new ArrayList<String>(attributes.size());
        for ( String key : header ) {
            if ( containsKey(key) ) { // avoid the header
                strings.add(this.get(key));
                //System.out.printf("Adding %s%n", this.get(key));
            }
        }
        return Utils.join("\t", strings);
    }

    // ----------------------------------------------------------------------
    //
    // map functions
    //
    // ----------------------------------------------------------------------
    public int size()                               { return attributes.size(); }
    public boolean isEmpty()                        { return attributes.isEmpty(); }
    public boolean containsValue(Object o)          { return attributes.containsValue(o); }
    public String remove(Object o)                  { return attributes.remove(o); }
    public void clear()                             { attributes.clear(); }
    public java.util.Set<String> keySet()           { return attributes.keySet(); }
    public java.util.Collection<String> values()    { return attributes.values(); }

    public void putAll(java.util.Map<? extends String, ? extends String> map) {
        attributes.putAll(map);
    }

    public java.util.Set<java.util.Map.Entry<String,String>> entrySet() {
        return attributes.entrySet();
    }

    // ----------------------------------------------------------------------
    //
    // formatting
    //
    // ----------------------------------------------------------------------
    public String toString() {
        return String.format("%s\t%s", getLocation(), getAttributeString());
    }

    /**
     * The delimiter regular expression that should be used to separate fields in data rows
     * and header.
     * 
     * @return
     */
    public String delimiterRegex() {
        return DELIMITER_REGEX;
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
        header = (ArrayList<String>)(headerObj);

        //System.out.printf("parts [len=%d] is '%s'%n", parts.length, Utils.join(":", parts));

        if ( parts.length == 0 || COMMENT_PATTERN.matcher(parts[0]).matches() || HEADER_PATTERN.matcher(parts[0]).matches() )
            return false;

        if (header.size() != parts.length) {
            throw new IOException(String.format("Header length %d not equal to Tabular parts length %d", header.size(), parts.length));
        }

        loc = GenomeLoc.parseGenomeLoc(parts[0]);
        for ( int i = 1; i < parts.length; i++ ) {
            put(header.get(i), parts[i]);
        }

        //System.out.printf("Parsed %s%n", this);

        return true;
    }
}