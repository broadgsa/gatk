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
 * To change this template use File | Settings | File Templates.
 */
public class TabularROD extends BasicReferenceOrderedDatum implements Map<String, String> {
    private static Logger logger = Logger.getLogger(TabularROD.class);

    private GenomeLoc loc;
    private HashMap<String, String> attributes;
    private ArrayList<String> header;

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
                header = new ArrayList<String>(Arrays.asList(line.split("\\s+")));
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

        //System.exit(1);
        return header;
    }

    private static int MAX_LINES_TO_LOOK_FOR_HEADER = 1000;
    private static Pattern HEADER_PATTERN = Pattern.compile("^\\s*HEADER.*");
    private static Pattern COMMENT_PATTERN = Pattern.compile("^#.*");

    // ----------------------------------------------------------------------
    //
    // Accessors
    //
    // ----------------------------------------------------------------------
    public GenomeLoc getLocation() {
        return loc;
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
                System.out.printf("Adding %s%n", this.get(key));
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

    public String delimiter() {
        return "\\s+";
    }

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