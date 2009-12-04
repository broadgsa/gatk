package org.broadinstitute.sting.utils;

import java.util.List;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * Common utilities for dealing with strings.
 *
 * @author mhanna
 * @version 0.1
 */
public class TextFormattingUtils {
    /**
     * The default line width, for GATK output written to the screen.
     */
    public static final int DEFAULT_LINE_WIDTH = 100;

    /**
     * Simple implementation of word-wrap for a line of text.  Idea and
     * regexp shamelessly stolen from http://joust.kano.net/weblog/archives/000060.html.
     * Regexp can probably be simplified for our application.
     * @param text Text to wrap.
     * @param width Maximum line width.
     * @return A list of word-wrapped lines.
     */
    public static List<String> wordWrap( String text, int width ) {
        Pattern wrapper = Pattern.compile( String.format(".{0,%d}(?:\\S(?: |$)|$)", width-1) );
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

}
