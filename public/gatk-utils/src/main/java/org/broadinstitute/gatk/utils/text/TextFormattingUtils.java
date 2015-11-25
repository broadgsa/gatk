/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.text;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.io.Resource;

import java.io.IOException;
import java.io.StringReader;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Common utilities for dealing with text formatting.
 *
 * @author mhanna
 * @version 0.1
 */
public class TextFormattingUtils {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(TextFormattingUtils.class);

    /**
     * The contents of the GATK bundle.  If no such resource exists, warn the user and create an empty bundle.
     */
    public static final ResourceBundle GATK_RESOURCE_BUNDLE = loadResourceBundle("GATKText", null);

    /**
     * The default line width, for GATK output written to the screen.
     */
    public static final int DEFAULT_LINE_WIDTH = 120;

    /**
     * Simple implementation of word-wrap for a line of text.  Idea and
     * regexp shamelessly stolen from http://joust.kano.net/weblog/archives/000060.html.
     * Regexp can probably be simplified for our application.
     * @param text Text to wrap.
     * @param width Maximum line width.
     * @return A list of word-wrapped lines.
     */
    public static List<String> wordWrap( String text, int width ) {
        Pattern wrapper = Pattern.compile( String.format(".{0,%d}(?:\\S(?:[\\s|]|$)|$)", width-1) );
        Matcher matcher = wrapper.matcher( text );

        List<String> wrapped = new ArrayList<String>();
        while( matcher.find() ) {
            // Regular expression is supersensitive to whitespace.
            // Assert that content is present before adding the line.
            String line = matcher.group().trim();
            if( line.length() > 0 )
                wrapped.add( matcher.group() );
        }
        return wrapped;
    }

    /**
     * Compares two strings independently of case sensitivity.
     */
    public static class CaseInsensitiveComparator implements Comparator<String> {
        /**
         * Compares the order of lhs to rhs, not taking case into account.
         * @param lhs First object to compare.
         * @param rhs Second object to compare.
         * @return 0 if objects are identical; -1 if lhs is before rhs, 1 if rhs is before lhs.  Nulls are treated as after everything else.
         */
        public int compare(String lhs, String rhs) {
            if(lhs == null && rhs == null) return 0;
            if(lhs == null) return 1;
            if(rhs == null) return -1;
            return lhs.toLowerCase().compareTo(rhs.toLowerCase());
        }
    }

    /**
     * Load the contents of a resource bundle with the given name.  If no such resource exists, warn the user
     * and create an empty bundle.
     * @param bundleName The name of the bundle to load.
     * @param relativeClass The relative class or null to load a bundle from the root.
     * @return The best resource bundle that can be found matching the given name.
     */
    public static ResourceBundle loadResourceBundle(String bundleName, Class<?> relativeClass) {
        final ResourceBundle.Control c = ResourceBundle.Control.getControl(ResourceBundle.Control.FORMAT_DEFAULT);
        final String resourceName = c.toResourceName(c.toBundleName(bundleName, Locale.ROOT), "properties");
        final Resource resource = new Resource(resourceName, relativeClass);
        ResourceBundle bundle;
        try {
            bundle = new PropertyResourceBundle(resource.getAllResourcesContentsAsStream());
        }
        catch(Exception ex) {
            //logger.warn("Unable to load help text.  Help output will be sparse.");
            // Generate an empty resource bundle.
            try {
                bundle = new PropertyResourceBundle(new StringReader(""));
            }
            catch(IOException ioe) {
                throw new ReviewedGATKException("No resource bundle found, and unable to create an empty placeholder.",ioe);
            }
        }
        return bundle;
    }


    /**
     * Returns the word starting positions within line, excluding the first position 0.
     * The returned list is compatible with splitFixedWidth.
     * @param line Text to parse.
     * @return the word starting positions within line, excluding the first position 0.
     */
    public static List<Integer> getWordStarts(String line) {
        if (line == null)
            throw new ReviewedGATKException("line is null");
        List<Integer> starts = new ArrayList<Integer>();
        int stop = line.length();
        for (int i = 1; i < stop; i++)
            if (Character.isWhitespace(line.charAt(i-1)))
                if(!Character.isWhitespace(line.charAt(i)))
                    starts.add(i);
        return starts;
    }

    /**
     * Parses a fixed width line of text.
     * @param line Text to parse.
     * @param columnStarts the column starting positions within line, excluding the first position 0.
     * @return The parsed string array with each entry trimmed.
     */
    public static String[] splitFixedWidth(String line, List<Integer> columnStarts) {
        if (line == null)
            throw new ReviewedGATKException("line is null");
        if (columnStarts == null)
            throw new ReviewedGATKException("columnStarts is null");
        int startCount = columnStarts.size();
        String[] row = new String[startCount + 1];
        if (startCount == 0) {
            row[0] = line.trim();
        } else {
            row[0] = line.substring(0, columnStarts.get(0)).trim();
            for (int i = 1; i < startCount; i++)
                row[i] = line.substring(columnStarts.get(i - 1), columnStarts.get(i)).trim();
            row[startCount] = line.substring(columnStarts.get(startCount - 1)).trim();
        }
        return row;
    }

    /**
     * Parses a line of text by whitespace.
     * @param line Text to parse.
     * @return The parsed string array.
     */
    public static String[] splitWhiteSpace(String line) {
        if (line == null)
            throw new ReviewedGATKException("line is null");
        return line.trim().split("\\s+");
    }
}
