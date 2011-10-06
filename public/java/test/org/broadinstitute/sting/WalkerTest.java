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
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.Tribble;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;

import java.io.File;
import java.util.*;

public class WalkerTest extends BaseTest {
    private static final boolean ENABLE_REPORTING = false;

    @BeforeMethod
    public void initializeRandomGenerator() {
        GenomeAnalysisEngine.resetRandomGenerator();
    }

    public MD5DB.MD5Match assertMatchingMD5(final String name, final File resultsFile, final String expectedMD5) {
        return MD5DB.assertMatchingMD5(name, resultsFile, expectedMD5, parameterize());
    }

    public void maybeValidateSupplementaryFile(final String name, final File resultFile) {
        File indexFile = Tribble.indexFile(resultFile);
        //System.out.println("Putative index file is " + indexFile);
        if ( indexFile.exists() ) {
            if ( resultFile.getAbsolutePath().contains(".vcf") ) {
                // todo -- currently we only understand VCF files! Blow up since we can't test them
                throw new StingException("Found an index created for file " + resultFile + " but we can only validate VCF files.  Extend this code!");
            }

            assertOnDiskIndexEqualToNewlyCreatedIndex(indexFile, name, resultFile);
        }
    }


    public static void assertOnDiskIndexEqualToNewlyCreatedIndex(final File indexFile, final String name, final File resultFile) {
        System.out.println("Verifying on-the-fly index " + indexFile + " for test " + name + " using file " + resultFile);
        Index indexFromOutputFile = IndexFactory.createIndex(resultFile, new VCFCodec());
        Index dynamicIndex = IndexFactory.loadIndex(indexFile.getAbsolutePath());

        if ( ! indexFromOutputFile.equalsIgnoreProperties(dynamicIndex) ) {
            Assert.fail(String.format("Index on disk from indexing on the fly not equal to the index created after the run completed.  FileIndex %s vs. on-the-fly %s%n",
                    indexFromOutputFile.getProperties(),
                    dynamicIndex.getProperties()));
        }
    }

    public List<String> assertMatchingMD5s(final String name, List<File> resultFiles, List<String> expectedMD5s) {
        List<String> md5s = new ArrayList<String>();
        List<MD5DB.MD5Match> fails = new ArrayList<MD5DB.MD5Match>();

        for (int i = 0; i < resultFiles.size(); i++) {
            MD5DB.MD5Match result = assertMatchingMD5(name, resultFiles.get(i), expectedMD5s.get(i));
            if ( ! result.failed ) {
                maybeValidateSupplementaryFile(name, resultFiles.get(i));
                md5s.add(result.md5);
            } else {
                fails.add(result);
            }
        }

        if ( ! fails.isEmpty() ) {
            for ( final MD5DB.MD5Match fail : fails ) {
                logger.warn("Fail: " + fail.failMessage);
            }
            Assert.fail("Test failed: " + name);
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
        String args = "";
        int nOutputFiles = -1;
        List<String> md5s = null;
        List<String> exts = null;
        Class expectedException = null;

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
        MD5DB.ensureMd5DbDirectory(); // ensure the md5 directory exists

        List<File> tmpFiles = new ArrayList<File>();
        for (int i = 0; i < spec.nOutputFiles; i++) {
            String ext = spec.exts == null ? ".tmp" : "." + spec.exts.get(i);
            File fl = createTempFile(String.format("walktest.tmp_param.%d", i), ext);
            tmpFiles.add(fl);
        }

        final String args = String.format(spec.args, tmpFiles.toArray());
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
    public static void executeTest(String name, String args, Class expectedException) {
        CommandLineGATK instance = new CommandLineGATK();
        String[] command = Utils.escapeExpressions(args);

        // add the logging level to each of the integration test commands
        command = Utils.appendArray(command, "-et", ENABLE_REPORTING ? "STANDARD" : "NO_ET");

        // run the executable
        boolean gotAnException = false;
        try {
            System.out.println(String.format("Executing test %s with GATK arguments: %s", name, Utils.join(" ",command)));
            CommandLineExecutable.start(instance, command);
        } catch (Exception e) {
            gotAnException = true;
            if ( expectedException != null ) {
                // we expect an exception
                System.out.println(String.format("Wanted exception %s, saw %s", expectedException, e.getClass()));
                if ( expectedException.isInstance(e) ) {
                    // it's the type we expected
                    System.out.println(String.format("  => %s PASSED", name));
                } else {
                    e.printStackTrace();
                    Assert.fail(String.format("Test %s expected exception %s but got %s instead",
                            name, expectedException, e.getClass()));
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
