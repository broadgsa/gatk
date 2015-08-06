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

package org.broadinstitute.gatk.engine.walkers;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.bcf2.BCF2Utils;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.gatk.engine.CommandLineExecutable;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.crypt.CryptUtils;
import org.broadinstitute.gatk.engine.phonehome.GATKRunReport;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.MD5DB;
import org.broadinstitute.gatk.utils.MD5Mismatch;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.classloader.JVMUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.testng.Assert;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.BeforeMethod;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.*;

public class WalkerTest extends BaseTest {
    public static final String gatkKeyFile = CryptUtils.GATK_USER_KEY_DIRECTORY + "gsamembers_broadinstitute.org.key";

    private static final boolean GENERATE_SHADOW_BCF = true;
    private static final boolean ENABLE_PHONE_HOME_FOR_TESTS = false;
    private static final boolean ENABLE_ON_THE_FLY_CHECK_FOR_VCF_INDEX = false;
    private static final boolean ENABLE_AUTO_INDEX_CREATION_AND_LOCKING_FOR_TESTS = false;

    private static MD5DB md5DB = new MD5DB();

    @BeforeMethod
    public void initializeWalkerTests() {
        logger.debug("Initializing walker tests");
        Utils.resetRandomGenerator();
    }

    @AfterSuite
    public void finalizeWalkerTests() {
        logger.debug("Finalizing walker tests");
        md5DB.close();
    }

    public static MD5DB getMd5DB() {
        return md5DB;
    }

    public void validateOutputBCFIfPossible(final String name, final File resultFile) {
        final File bcfFile = BCF2Utils.shadowBCF(resultFile);
        if ( bcfFile != null && bcfFile.exists() ) {
            logger.warn("Checking shadow BCF output file " + bcfFile + " against VCF file " + resultFile);
            try {
                assertVCFandBCFFilesAreTheSame(resultFile, bcfFile);
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
                throw new GATKException("Found an index created for file " + resultFile + " but we can only validate VCF files.  Extend this code!");
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

    public List<String> assertMatchingMD5s(final String testName, final String testClassName, List<File> resultFiles, List<String> expectedMD5s) {
        List<String> md5s = new ArrayList<String>();
        List<MD5DB.MD5Match> fails = new ArrayList<MD5DB.MD5Match>();

        for (int i = 0; i < resultFiles.size(); i++) {
            MD5DB.MD5Match result = getMd5DB().testFileMD5(testName, testClassName, resultFiles.get(i), expectedMD5s.get(i), parameterize());
            validateOutputBCFIfPossible(testName, resultFiles.get(i));
            if ( ! result.failed ) {
                validateOutputIndex(testName, resultFiles.get(i));
                md5s.add(result.expectedMD5);
            } else {
                fails.add(result);
            }
        }

        if ( ! fails.isEmpty() ) {
            List<String> actuals = new ArrayList<String>();
            List<String> expecteds = new ArrayList<String>();
            List<String> diffEngineOutputs = new ArrayList<String>();

            for ( final MD5DB.MD5Match fail : fails ) {
                actuals.add(fail.actualMD5);
                expecteds.add(fail.expectedMD5);
                diffEngineOutputs.add(fail.diffEngineOutput);
                logger.warn("Fail: " + fail.failMessage);
            }

            final MD5Mismatch failure = new MD5Mismatch(actuals, expecteds, diffEngineOutputs);
            Assert.fail(failure.toString());
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

        // Name of the test class that created this test case
        private Class testClass;

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
            this.testClass = getCallingTestClass();
        }

        public WalkerTestSpec(String args, List<String> exts, List<String> md5s) {
            this(args, -1, exts, md5s);
        }

        public WalkerTestSpec(String args, int nOutputFiles, List<String> exts, List<String> md5s) {
            this.args = args;
            this.nOutputFiles = md5s.size();
            this.md5s = md5s;
            this.exts = exts;
            this.testClass = getCallingTestClass();
        }

        // @Test(expectedExceptions) doesn't work in integration tests, so use this instead
        public WalkerTestSpec(String args, int nOutputFiles, Class expectedException) {
            this.args = args;
            this.nOutputFiles = nOutputFiles;
            this.expectedException = expectedException;
            this.testClass = getCallingTestClass();
        }

        private Class getCallingTestClass() {
            return JVMUtils.getCallingClass(getClass());
        }

        public String getTestClassName() {
            return testClass.getSimpleName();
        }

        public String getArgsWithImplicitArgs() {
            String args = this.args;
            if ( includeImplicitArgs ) {
                args = args + (ENABLE_PHONE_HOME_FOR_TESTS ?
                        String.format(" -et %s ", GATKRunReport.PhoneHomeOption.AWS) :
                        String.format(" -et %s -K %s ", GATKRunReport.PhoneHomeOption.NO_ET, gatkKeyFile));
                if ( includeShadowBCF && GENERATE_SHADOW_BCF )
                    args = args + " --generateShadowBCF ";
                if ( ! ENABLE_AUTO_INDEX_CREATION_AND_LOCKING_FOR_TESTS )
                    args = args + " --disable_auto_index_creation_and_locking_when_reading_rods ";
            }

            return args;
        }

        /**
         * In the case where the input VCF files are malformed and cannot be fixed
         * this function tells the engine to not try to generate a shadow BCF
         * which will ultimately blow up...
         */
        public void disableShadowBCF() { this.includeShadowBCF = false; }
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
            if ( ! expectsException() ) throw new ReviewedGATKException("Tried to get expection for walker test that doesn't expect one");
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

    public enum ParallelTestType {
        TREE_REDUCIBLE,
        NANO_SCHEDULED,
        BOTH
    }

    protected Pair<List<File>, List<String>> executeTestParallel(final String name, WalkerTestSpec spec, ParallelTestType testType) {
        final List<Integer> ntThreads  = testType == ParallelTestType.TREE_REDUCIBLE || testType == ParallelTestType.BOTH ? Arrays.asList(1, 4) : Collections.<Integer>emptyList();
        final List<Integer> cntThreads = testType == ParallelTestType.NANO_SCHEDULED || testType == ParallelTestType.BOTH ? Arrays.asList(1, 4) : Collections.<Integer>emptyList();

        return executeTest(name, spec, ntThreads, cntThreads);
    }

    protected Pair<List<File>, List<String>> executeTestParallel(final String name, WalkerTestSpec spec) {
        return executeTestParallel(name, spec, ParallelTestType.TREE_REDUCIBLE);
    }

    protected Pair<List<File>, List<String>> executeTest(final String name, WalkerTestSpec spec, List<Integer> ntThreads, List<Integer> cpuThreads) {
        String originalArgs = spec.args;
        Pair<List<File>, List<String>> results = null;

        boolean ran1 = false;
        for ( int nt : ntThreads ) {
            String extra = nt == 1 ? "" : (" -nt " + nt);
            ran1 = ran1 || nt == 1;
            spec.args = originalArgs + extra;
            results = executeTest(name + "-nt-" + nt, spec);
        }

        for ( int nct : cpuThreads ) {
            if ( nct != 1 ) {
                String extra = " -nct " + nct;
                spec.args = originalArgs + extra;
                results = executeTest(name + "-cnt-" + nct, spec);
            }
        }

        return results;
    }

    protected Pair<List<File>, List<String>> executeTest(final String name, WalkerTestSpec spec) {
        List<File> tmpFiles = new ArrayList<File>();
        for (int i = 0; i < spec.nOutputFiles; i++) {
            String ext = spec.exts == null ? ".tmp" : "." + spec.exts.get(i);
            File fl = createTempFile(String.format("walktest.tmp_param.%d", i), ext);

            // Cleanup any potential shadow BCFs on exit too, if we're generating them
            if ( spec.includeShadowBCF && GENERATE_SHADOW_BCF ) {
                final File potentalShadowBCFFile = BCF2Utils.shadowBCF(fl);
                potentalShadowBCFFile.deleteOnExit();
                new File(potentalShadowBCFFile.getAbsolutePath() + Tribble.STANDARD_INDEX_EXTENSION).deleteOnExit();
            }

            tmpFiles.add(fl);
        }

        final String args = String.format(spec.getArgsWithImplicitArgs(), tmpFiles.toArray());
        System.out.println(Utils.dupString('-', 80));

        if ( spec.expectsException() ) {
            // this branch handles the case were we are testing that a walker will fail as expected
            return executeTest(name, spec.getTestClassName(), spec.getOutputFileLocation(), null, tmpFiles, args, spec.getExpectedException());
        } else {
            List<String> md5s = new LinkedList<String>();
            md5s.addAll(spec.md5s);

            // check to see if they included any auxillary files, if so add them to the list and set them to be deleted on exit
            for (String md5 : spec.auxillaryFiles.keySet()) {
                md5s.add(md5);
                final File auxFile = spec.auxillaryFiles.get(md5);
                auxFile.deleteOnExit();
                tmpFiles.add(auxFile);
            }
            return executeTest(name, spec.getTestClassName(), spec.getOutputFileLocation(), md5s, tmpFiles, args, null);
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
     * @param testName     the name of the test
     * @param testClassName the name of the class that contains the test
     * @param md5s     the list of md5s
     * @param tmpFiles the temp file corresponding to the md5 list
     * @param args     the argument list
     * @param expectedException the expected exception or null
     * @return a pair of file and string lists
     */
    private Pair<List<File>, List<String>> executeTest(String testName, String testClassName, File outputFileLocation, List<String> md5s, List<File> tmpFiles, String args, Class expectedException) {
        if ( md5s != null ) qcMD5s(testName, md5s);

        if (outputFileLocation != null)
            args += " -o " + outputFileLocation.getAbsolutePath();
        executeTest(testName, testClassName, args, expectedException);

        if ( expectedException != null ) {
            return null;
        } else {
            // we need to check MD5s
            return new Pair<List<File>, List<String>>(tmpFiles, assertMatchingMD5s(testName, testClassName, tmpFiles, md5s));
        }
    }
    
    /**
     * execute the test, given the following:
     * @param testName      the name of the test
     * @param testClassName the name of the class that contains the test
     * @param args          the argument list
     * @param expectedException the expected exception or null
     */
    private void executeTest(String testName, String testClassName, String args, Class expectedException) {
        CommandLineGATK instance = new CommandLineGATK();
        String[] command = Utils.escapeExpressions(args);
        // run the executable
        boolean gotAnException = false;
        try {
            final String now = new SimpleDateFormat("HH:mm:ss").format(new Date());
            final String cmdline = Utils.join(" ",command);
            System.out.println(String.format("[%s] Executing test %s:%s with GATK arguments: %s", now, testClassName, testName, cmdline));
            // also write the command line to the HTML log for convenient follow-up
            // do the replaceAll so paths become relative to the current
            BaseTest.log(cmdline.replaceAll(publicTestDirRoot, "").replaceAll(privateTestDirRoot, ""));
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
                    final String message = String.format("Test %s:%s expected exception %s but instead got %s with error message %s",
                            testClassName, testName, expectedException, e.getClass(), e.getMessage());
                    if ( e.getCause() != null ) {
                        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
                        final PrintStream ps = new PrintStream(baos);
                        e.getCause().printStackTrace(ps);
                        BaseTest.log(message);
                        BaseTest.log(baos.toString());
                    }
                    Assert.fail(message);
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
                Assert.fail(String.format("Test %s:%s expected exception %s but none was thrown", testClassName, testName, expectedException.toString()));
        } else {
            if ( CommandLineExecutable.result != 0) {
                throw new RuntimeException("Error running the GATK with arguments: " + args);
            }
        }
    }


    protected File createTempFileFromBase(final String name) {
        File fl = new File(name);
        fl.deleteOnExit();
        return fl;
    }
}
