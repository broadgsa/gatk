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

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ListAppenderUnitTest {

    private final static int EVENTS_TO_KEEP = 10;
    private final static int EVENTS = 200;
    private final static String MSG_TO_LOG = ": All work and no play makes Jack a dull boy.";
    private final static String MSG_TOO_LONG = ": Don't write this one.";
    private final static Logger LOGGER = Logger.getLogger(ListAppenderUnitTest.class);

    @Test
    public void testListAppender() {

        LOGGER.removeAllAppenders();
        final ListAppender listAppender = new ListAppender(new PatternLayout(PatternLayout.DEFAULT_CONVERSION_PATTERN), EVENTS_TO_KEEP, Level.WARN);
        LOGGER.addAppender(listAppender);

        for ( int i = 0; i < EVENTS; i++ ) {
            LOGGER.warn(i + MSG_TO_LOG);
            LOGGER.info(i + MSG_TOO_LONG);
        }

        final String listAppenderString = listAppender.toString();
        listAppender.write();

        Assert.assertEquals(StringUtils.countMatches(listAppenderString, MSG_TO_LOG), EVENTS_TO_KEEP);
        Assert.assertFalse(listAppenderString.contains(MSG_TOO_LONG));
        Assert.assertTrue(listAppenderString.contains(Integer.toString(EVENTS_TO_KEEP)));
        Assert.assertTrue(listAppenderString.contains(listAppender.logLevelToKeep.toString()));
        Assert.assertEquals(listAppender.numLoggingEvents, EVENTS);
    }

    @Test
    public void testListAppenderNoWarnMsgs() {

        LOGGER.removeAllAppenders();
        final ListAppender listAppender = new ListAppender(new PatternLayout(PatternLayout.DEFAULT_CONVERSION_PATTERN), EVENTS_TO_KEEP, Level.WARN);
        LOGGER.addAppender(listAppender);

        for ( int i = 0; i < EVENTS; i++ ) {
            LOGGER.info(i + MSG_TOO_LONG);
        }

        final String listAppenderString = listAppender.toString();
        Assert.assertEquals(StringUtils.countMatches(listAppenderString, MSG_TO_LOG), 0);
        Assert.assertFalse(listAppenderString.contains(MSG_TOO_LONG));
        Assert.assertFalse(listAppenderString.contains(Integer.toString(EVENTS_TO_KEEP)));
        Assert.assertTrue(listAppenderString.contains(listAppender.logLevelToKeep.toString().toLowerCase()));
        Assert.assertEquals(listAppender.numLoggingEvents, 0);
    }

    @Test(expectedExceptions = ReviewedGATKException.class)
    public void testListAppenderMissingLayout() {

        LOGGER.removeAllAppenders();
        final ListAppender listAppender = new ListAppender(null, 1, Level.WARN);
        LOGGER.addAppender(listAppender);
        LOGGER.warn("Warning");

        final String listAppenderString = listAppender.toString();
    }
}
