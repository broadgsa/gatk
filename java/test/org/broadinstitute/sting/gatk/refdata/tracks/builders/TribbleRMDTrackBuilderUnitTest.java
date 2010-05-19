/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.tracks.builders;


import org.apache.log4j.AppenderSkeleton;
import org.apache.log4j.Logger;
import org.apache.log4j.spi.LoggingEvent;
import org.broad.tribble.vcf.VCFCodec;
import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class TribbleRMDTrackBuilderUnitTest
 *
 * Testing out the builder for tribble Tracks
 */
public class TribbleRMDTrackBuilderUnitTest extends BaseTest {
    private TribbleRMDTrackBuilder builder;

    @Before
    public void setup() {
        builder = new TribbleRMDTrackBuilder();
    }

    @Test
    public void testBuilder() {
        Map<String,Class> classes = builder.getAvailableTrackNamesAndTypes();
        Assert.assertTrue(classes.size() > 0);
    }

    @Test
    public void testBuilderIndexUnwriteable() {
        File vcfFile = new File(validationDataLocation + "/ROD_validation/relic.vcf");
        try {
            builder.loadIndex(vcfFile,new VCFCodec(), true);
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        // make sure we didn't write the file (check that it's timestamp is within bounds)
        Assert.assertTrue(Math.abs(1274210993000l - new File(vcfFile + TribbleRMDTrackBuilder.linearIndexExtension).lastModified()) < 100);

    }

    @Test
    public void testBuilderIndexOutOfDate() {
        Logger logger = Logger.getLogger(TribbleRMDTrackBuilder.class);
        ValidationAppender appender = new ValidationAppender("Tribble index file /humgen/gsa-hpprojects/GATK/data/Validation_Data/ROD_validation/newerTribbleTrack.vcf.idx is older than the track file /humgen/gsa-hpprojects/GATK/data/Validation_Data/ROD_validation/newerTribbleTrack.vcf, this can lead to unexpected behavior");
        logger.addAppender(appender);
        File vcfFile = new File(validationDataLocation + "/ROD_validation/newerTribbleTrack.vcf");
        try {
            builder.loadIndex(vcfFile,new VCFCodec(), true);
        } catch (IOException e) {
            e.printStackTrace();
            Assert.fail("IO exception unexpected" + e.getMessage());
        }
        // check to make sure the appender saw the target string 
        Assert.assertTrue(appender.foundString());

    }
}

/**
 * this appender looks for a specific message in the log4j stream.
 * It can be used to verify that a specific message was generated to the logging system.
 */
class ValidationAppender extends AppenderSkeleton {

    private boolean foundString = false;
    private String targetString = "";

    public ValidationAppender(String target) {
        targetString = target;
    }

    @Override
    protected void append(LoggingEvent loggingEvent) {
        if (loggingEvent.getMessage().equals(targetString))
            foundString = true;
    }

    @Override
    public void close() {
        // do nothing
    }

    @Override
    public boolean requiresLayout() {
        return false;
    }

    public boolean foundString() {
        return foundString;
    }
}
