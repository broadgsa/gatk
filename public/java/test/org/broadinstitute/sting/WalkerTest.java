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

package org.broadinstitute.sting;

import org.apache.commons.lang.StringUtils;
import org.broad.tribble.Tribble;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContextTestProvider;

import java.io.*;

import org.testng.Assert;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.BeforeMethod;

import java.text.SimpleDateFormat;
import java.util.*;

public class WalkerTest extends BaseTest {
    private static final boolean GENERATE_SHADOW_BCF = false;
    private static final boolean ENABLE_PHONE_HOME_FOR_TESTS = false;
    private static final boolean ENABLE_ON_THE_FLY_CHECK_FOR_VCF_INDEX = false;

    private static MD5DB md5DB = new MD5DB();

    @BeforeMethod
    public void initializeWalkerTests() {
        logger.debug("Initializing walker tests");
        GenomeAnalysisEngine.resetRandomGenerator();
    }

    @AfterSuite
    public void finalizeWalkerTests() {
        logger.debug("Finalizing walker tests");
        md5DB.close();
    }

    public static MD5DB getMd5DB() {
        return md5DB;
    }

    public MD5DB.MD5Match assertMatchingMD5(final String name, final File resultsFile, final String expectedMD5) {
        return getMd5DB().assertMatchingMD5(name, resultsFile, expectedMD5, parameterize());
    }

    public void validateOutputBCFIfPossible(final String name, final File resultFile) {
        final File bcfFile = BCF2Utils.shadowBCF(resultFile);
        if ( bcfFile.exists() ) {
            logger.warn("Checking shadow BCF output file " + bcfFile + " against VCF file " + resultFile);
            try {
                VariantContextTestProvider.assertVCFandBCFFilesAreTheSame(resultFile, bcfFile);
                logger.warn("  Shadow BCF PASSED!");
            } catch ( Exception e ) {
                Assert.fail("Exception received reading shadow BCFFile " + bcfFile + " for test " + name, e);
            }
        }
    }

    public void validateOutputIndex(final String name, final File resultFile) {
        if ( !ENABLE_ON_THE_FLY_CHECK_FOR_VCF_INDEX )
            return;

        File indexFile = Tribble.indexFile(resultFile);
        //System.out.println("Putative index file is " + indexFile);
        if ( indexFile.exists() ) {
            if ( resultFile.getAbsolutePath().contains(".vcf") ) {
                // todo -- currently we only understand VCF files! Blow up since we can't test them
                throw new StingException("Found an index created for file " + resultFile + " but we can only validate VCF files.  Extend this code!");
            }

            System.out.println("Verifying on-the-fly index " + indexFile + " for test " + name + " using file " + resultFile);
            Index indexFromOutputFile = IndexFactory.createDynamicIndex(resultFile, new VCFCodec());
            Index dynamicIndex = IndexFactory.loadIndex(indexFile.getAbsolutePath());

            if ( ! indexFromOutputFile.equalsIgnoreProperties(dynamicIndex) ) {
                Assert.fail(String.format("Index on disk from indexing on the fly not equal to the index created after the run completed.  FileIndex %s vs. on-the-fly %s%n",
                        indexFromOutputFile.getProperties(),
                        dynamicIndex.getProperties()));
            }
        }
    }

    public List<String> assertMatchingMD5s(final String name, List<File> resultFiles, List<String> expectedMD5s) {
        List<String> md5s = new ArrayList<String>();
        List<MD5DB.MD5Match> fails = new ArrayList<MD5DB.MD5Match>();

        for (int i = 0; i < resultFiles.size(); i++) {
            MD5DB.MD5Match result = assertMatchingMD5(name, resultFiles.get(i), expectedMD5s.get(i));
            validateOutputBCFIfPossible(name, resultFiles.get(i));
            if ( ! result.failed ) {
                validateOutputIndex(name, resultFiles.get(i));
                md5s.add(result.expectedMD5);
            } else {
                fails.add(result);
            }
        }

        if ( ! fails.isEmpty() ) {
            List<String> actuals = new ArrayList<String>();
            List<String> expecteds = new ArrayList<String>();
            for ( final MD5DB.MD5Match fail : fails ) {
                actuals.add(fail.actualMD5);
                expecteds.add(fail.expectedMD5);
                logger.warn("Fail: " + fail.failMessage);
            }

            final MD5Mismatch failure = new MD5Mismatch(actuals, expecteds);
            Assert.fail(failure.toString(), failure);
        }

        return md5s;
    }

    public String buildCommandLine(String... arguments) {
        String cmdline = "";

        for ( int argIndex = 0; argIndex < arguments.length; argIndex++ ) {
            cmdline += arguments[argIndex];

            if (argIndex < arguments.length - 1) {
                cmdline += " ";
            }
        }

        return cmdline;
    }

    public class WalkerTestSpec {
        // Arguments implicitly included in all Walker command lines, unless explicitly
        // disabled using the disableImplicitArgs() method below.
        String args = "";
        int nOutputFiles = -1;
        List<String> md5s = null;
        List<String> exts = null;
        Class expectedException = null;
        boolean includeImplicitArgs = true;
        boolean includeShadowBCF = true;
        boolean repairHeader = false;

        // the default output path for the integration test
        private File outputFileLocation = null;

        protected Map<String, File> auxillaryFiles = new HashMap<String, File>();

        public WalkerTestSpec(String args, List<String> md5s) {
            this(args, -1, md5s);
        }

        public WalkerTestSpec(String args, int nOutputFiles, List<String> md5s) {
            this.args = args;
            this.nOutputFiles = md5s.size();
            this.md5s = md5s;
        }

        public WalkerTestSpec(String args, List<String> exts, List<String> md5s) {
            this(args, -1, exts, md5s);
        }

        public WalkerTestSpec(String args, int nOutputFiles, List<String> exts, List<String> md5s) {
            this.args = args;
            this.nOutputFiles = md5s.size();
            this.md5s = md5s;
            this.exts = exts;
        }

        public WalkerTestSpec(String args, int nOutputFiles, Class expectedException) {
            this.args = args;
            this.nOutputFiles = nOutputFiles;
            this.expectedException = expectedException;
        }

        public String getArgsWithImplicitArgs() {
            String args = this.args;
            if ( includeImplicitArgs ) {
                args = args + (ENABLE_PHONE_HOME_FOR_TESTS ?
                        String.format(" -et %s ", GATKRunReport.PhoneHomeOption.STANDARD) :
                        String.format(" -et %s -K %s ", GATKRunReport.PhoneHomeOption.NO_ET, gatkKeyFile));
                if ( includeShadowBCF && GENERATE_SHADOW_BCF )
                    args = args + " --generateShadowBCF ";
                if ( repairHeader )
                    args = args + " --repairVCFHeader public/data/vcfHeaderForRepairs.vcf ";
            }

            return args;
        }

        /**
         * In the case where the input VCF files are malformed and cannot be fixed
         * this function tells the engine to not try to generate a shadow BCF
         * which will ultimately blow up...
         */
        public void disableShadowBCF() { this.includeShadowBCF = false; }
        public void repairHeaders() { this.repairHeader = true; }
        public void setOutputFileLocation(File outputFileLocation) {
            this.outputFileLocation = outputFileLocation;
        }        

        protected File getOutputFileLocation() {
            return outputFileLocation;
        }
        
        public boolean expectsException() {
            return expectedException != null;
        }

        public Class getExpectedException() {
            if ( ! expectsException() ) throw new ReviewedStingException("Tried to get expection for walker test that doesn't expect one");
            return expectedException;
        }

        public void addAuxFile(String expectededMD5sum, File outputfile) {
            auxillaryFiles.put(expectededMD5sum, outputfile);
        }

        public void disableImplicitArgs() {
            includeImplicitArgs = false;
        }
    }

    protected boolean parameterize() {
        return false;
    }

    protected Pair<List<File>, List<String>> executeTestParallel(final String name, WalkerTestSpec spec) {
        return executeTest(name, spec, Arrays.asList(1, 4));
    }

    protected Pair<List<File>, List<String>> executeTest(final String name, WalkerTestSpec spec, List<Integer> parallelThreads) {
        String originalArgs = spec.args;
        Pair<List<File>, List<String>> results = null;

        for ( int nt : parallelThreads ) {
            String extra = nt == 1 ? "" : (" -nt " + nt);
            spec.args = originalArgs + extra;
            results = executeTest(name + "-nt-" + nt, spec);
        }

        return results;
    }

    protected Pair<List<File>, List<String>> executeTest(final String name, WalkerTestSpec spec) {
        List<File> tmpFiles = new ArrayList<File>();
        for (int i = 0; i < spec.nOutputFiles; i++) {
            String ext = spec.exts == null ? ".tmp" : "." + spec.exts.get(i);
            File fl = createTempFile(String.format("walktest.tmp_param.%d", i), ext);
            tmpFiles.add(fl);
        }

        final String args = String.format(spec.getArgsWithImplicitArgs(), tmpFiles.toArray());
        System.out.println(Utils.dupString('-', 80));

        if ( spec.expectsException() ) {
            // this branch handles the case were we are testing that a walker will fail as expected
            return executeTest(name, spec.getOutputFileLocation(), null, tmpFiles, args, spec.getExpectedException());
        } else {
            List<String> md5s = new LinkedList<String>();
            md5s.addAll(spec.md5s);

            // check to see if they included any auxillary files, if so add them to the list
            for (String md5 : spec.auxillaryFiles.keySet()) {
                md5s.add(md5);
                tmpFiles.add(spec.auxillaryFiles.get(md5));
            }
            return executeTest(name, spec.getOutputFileLocation(), md5s, tmpFiles, args, null);
        }
    }

    private void qcMD5s(String name, List<String> md5s) {
        final String exampleMD5 = "709a1f482cce68992c637da3cff824a8";
        for (String md5 : md5s) {
            if ( md5 == null )
                throw new IllegalArgumentException("Null MD5 found in test " + name);
            if ( md5.equals("") ) // ok
                continue;
            if ( ! StringUtils.isAlphanumeric(md5) )
                throw new IllegalArgumentException("MD5 contains non-alphanumeric characters test " + name + " md5=" + md5);
            if ( md5.length() != exampleMD5.length() )
                throw new IllegalArgumentException("Non-empty MD5 of unexpected number of characters test " + name + " md5=" + md5);
        }
    }


    /**
     * execute the test, given the following:
     * @param name     the name of the test
     * @param md5s     the list of md5s
     * @param tmpFiles the temp file corresponding to the md5 list
     * @param args     the argument list
     * @param expectedException the expected exception or null
     * @return a pair of file and string lists
     */
    private Pair<List<File>, List<String>> executeTest(String name, File outputFileLocation, List<String> md5s, List<File> tmpFiles, String args, Class expectedException) {
        if ( md5s != null ) qcMD5s(name, md5s);

        if (outputFileLocation != null)
            args += " -o " + outputFileLocation.getAbsolutePath();
        executeTest(name, args, expectedException);

        if ( expectedException != null ) {
            return null;
        } else {
            // we need to check MD5s
            return new Pair<List<File>, List<String>>(tmpFiles, assertMatchingMD5s(name, tmpFiles, md5s));
        }
    }
    
    /**
     * execute the test, given the following:
     * @param name     the name of the test
     * @param args     the argument list
     * @param expectedException the expected exception or null
     */
    private void executeTest(String name, String args, Class expectedException) {
        CommandLineGATK instance = new CommandLineGATK();
        String[] command = Utils.escapeExpressions(args);

        // run the executable
        boolean gotAnException = false;
        try {
            final String now = new SimpleDateFormat("HH:mm:ss").format(new Date());
            final String cmdline = Utils.join(" ",command);
            System.out.println(String.format("[%s] Executing test %s with GATK arguments: %s", now, name, cmdline));
            // also write the command line to the HTML log for convenient follow-up
            // do the replaceAll so paths become relative to the current
            BaseTest.log(cmdline.replaceAll(testDirRoot, ""));
            CommandLineExecutable.start(instance, command);
        } catch (Exception e) {
            gotAnException = true;
            if ( expectedException != null ) {
                // we expect an exception
                //System.out.println(String.format("Wanted exception %s, saw %s", expectedException, e.getClass()));
                if ( expectedException.isInstance(e) ) {
                    // it's the type we expected
                    //System.out.println(String.format("  => %s PASSED", name));
                } else {
                    if ( e.getCause() != null )
                        e.getCause().printStackTrace(System.out);  // must print to stdout to see the message
                    Assert.fail(String.format("Test %s expected exception %s but instead got %s with error message %s",
                            name, expectedException, e.getClass(), e.getMessage()));
                }
            } else {
                // we didn't expect an exception but we got one :-(
                throw new RuntimeException(e);
            }
        }

        // catch failures from the integration test
        if ( expectedException != null ) {
            if ( ! gotAnException )
                // we expected an exception but didn't see it
                Assert.fail(String.format("Test %s expected exception %s but none was thrown", name, expectedException.toString()));
        } else {
            if ( CommandLineExecutable.result != 0) {
                throw new RuntimeException("Error running the GATK with arguments: " + args);
            }
        }
    }


    protected File createTempFileFromBase(String name) {
        File fl = new File(name);
        fl.deleteOnExit();
        return fl;
    }
}
