/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.utils.text;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.io.StringReader;
import java.io.IOException;

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
     * @return The best resource bundle that can be found matching the given name.
     */
    public static ResourceBundle loadResourceBundle(String bundleName) {
        ResourceBundle bundle;
        try {
            bundle = ResourceBundle.getBundle(bundleName);
        }
        catch(MissingResourceException ex) {
            logger.warn("Unable to load help text.  Help output will be sparse.");
            // Generate an empty resource bundle.
            try {
                bundle = new PropertyResourceBundle(new StringReader(""));
            }
            catch(IOException ioe) {
                throw new GATKException("No resource bundle found, and unable to create an empty placeholder.",ioe);
            }
        }
        return bundle;
    }

}
