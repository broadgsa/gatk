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

package org.broadinstitute.gatk.engine.arguments;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Level;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.MD5DB;
import org.broadinstitute.gatk.utils.MD5Mismatch;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.runtime.*;

public class LoggingIntegrationTest {
    private final MD5DB md5db = new MD5DB();

    private class LoggingTestProvider extends BaseTest.TestDataProvider {

        private final String baseCmdLine;

        private final Level logLevel;
        private final String logFileStr;
        public final File argumentOutputFile;
        public final File pipedOutputFile;

        private LoggingTestProvider(final Level logLevel, final boolean explicitLogfile) throws IOException {
            super(LoggingTestProvider.class);

            // TODO: a better command line that exercises log levels besides INFO
            this.baseCmdLine = String.format("java -cp %s %s -T TestPrintVariantsWalker -R %s -V %s -L 1:1000000-2000000 --no_cmdline_in_header",
                    StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                    CommandLineGATK.class.getCanonicalName(), BaseTest.b37KGReference, BaseTest.b37_NA12878_OMNI);

            this.logLevel = logLevel;
            this.logFileStr = explicitLogfile ? " -log " + BaseTest.createTempFile(logLevel.toString(), "log") : "";
            this.argumentOutputFile = BaseTest.createTempFile(logLevel.toString(), "vcf");
            this.pipedOutputFile = BaseTest.createTempFile(logLevel.toString(), "vcf");
        }

        public final String getCmdLine(boolean redirectStdout) {
            String command = String.format("%s -l %s %s", baseCmdLine, logLevel, logFileStr);
            return redirectStdout ? command : command + " -o " + argumentOutputFile;
        }

        public String toString() {
            return String.format("LoggingTestProvider logLevel=%s", logLevel);
        }
    }

    @DataProvider(name = "LoggingTest")
    public Object[][] makeLoggingTestProvider() throws IOException {
        for (Boolean explicitLogFile : Arrays.asList(true, false)) {
            // TODO: enable other logging levels when tests for those exist
            new LoggingTestProvider(Level.DEBUG, explicitLogFile);
        }

        return LoggingTestProvider.getTests(LoggingTestProvider.class);
    }

    /**
     * test that using an output argument produces the same output as stdout
     */
    @Test(dataProvider = "LoggingTest")
    public void testStdoutEquivalence(final LoggingTestProvider cfg) throws IOException {

        ProcessController pc = ProcessController.getThreadLocal();

        // output argument

        ProcessSettings ps = new ProcessSettings(cfg.getCmdLine(false).split("\\s+"));
        pc.execAndCheck(ps);
        String output_argument_md5 = md5db.calculateFileMD5(cfg.argumentOutputFile);

        // pipe to stdout

        ps = new ProcessSettings(cfg.getCmdLine(true).split("\\s+"));
        ps.setStdoutSettings(new OutputStreamSettings(cfg.pipedOutputFile));
        pc.execAndCheck(ps);

        MD5DB.MD5Match result = md5db.testFileMD5("LoggingIntegrationTest", "LoggingIntegrationTest", cfg.pipedOutputFile, output_argument_md5, false);
        if(result.failed) {
            final MD5Mismatch failure = new MD5Mismatch(result.actualMD5, result.expectedMD5, result.diffEngineOutput);
            Assert.fail(failure.toString());
        }
    }
}