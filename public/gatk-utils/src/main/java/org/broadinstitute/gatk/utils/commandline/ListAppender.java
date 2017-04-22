/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.commandline;

import org.apache.log4j.AppenderSkeleton;
import org.apache.log4j.Layout;
import org.apache.log4j.Level;
import org.apache.log4j.spi.LoggingEvent;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.ArrayList;
import java.util.List;

/**
 * Appender for saving logging loggingEvents to a list
 */
public class ListAppender extends AppenderSkeleton {
    /** 
     * List of logging events 
     *
     */
    protected final List<LoggingEvent> loggingEvents = new ArrayList<>();
    
    /** 
     * number of logging events to keep 
     */
    protected int loggingEventsToKeep = 0;

    /** 
     * number of logging events to keep 
     */
    protected int numLoggingEvents = 0;

    /**
     * Log warning level to keep
     */
    protected Level logLevelToKeep = Level.WARN;

    /**
     * Constructor
     *
     * @param layout               output layout
     * @param loggingEventsToKeep  number of warnings to keep
     * @param logLevelToKeep       logging logLevelToKeep ( TRACE, DEBUG, INFO, WARN, ERROR, FATAL, OFF )
     */
    public ListAppender(final Layout layout, final int loggingEventsToKeep, final Level logLevelToKeep) {
        this.layout = layout;
        this.loggingEventsToKeep = loggingEventsToKeep;
        this.logLevelToKeep = logLevelToKeep;
    }

    /**
     * Constructor
     * 
     * Hide from use
     */
    private ListAppender() {}

    /**
     * Close the appender
     *
     * Clear out the logging events
     */
    @Override
    public synchronized void close() {
        loggingEvents.clear();
        numLoggingEvents = 0;
    }

    /**
     * Is a message layout required?
     * 
     * @return true
     */
    @Override
    public boolean requiresLayout() {
        return true;
    }

    /**
     * Append log event to the list of logging events
     *
     * @param loggingEvent  The logging events
     */
    @Override
    protected synchronized void append(final LoggingEvent loggingEvent) {
        if ( loggingEvent.getLevel().equals(logLevelToKeep) ) {
            numLoggingEvents++;
            if ( numLoggingEvents <= loggingEventsToKeep )
                loggingEvents.add(loggingEvent);
        }
    }

    /**
     * The string representation for the class data
     * 
     * @return string containing the logging events
     * @throws ReviewedGATKException if layout is null
     */
    public synchronized String toString() {
        final String msgType = logLevelToKeep.toString() + " messages";
        if ( loggingEvents.isEmpty() ) {
            return "There were no " + new String(msgType).toLowerCase() + ".\n";
        } else {
            String out = "There were " + Integer.toString(numLoggingEvents) + " " + msgType;
            if ( layout == null )
                throw new ReviewedGATKException("layout cannot be null");
            out += ", the first " + loggingEvents.size() + " are repeated below.\n";
            for ( LoggingEvent event : loggingEvents ) {
                out += layout.format(event);
            }
            return out;
        }
    }

    /**
     * Write class data to standard out
     */
    public void write() {
        System.out.print(this);
    }

}
