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

import junit.framework.Assert;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.Utils;
import org.junit.Test;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.util.*;

public class WalkerTest extends BaseTest {
    // the default output path for the integration test
    private File outputFileLocation = null;
    private static final boolean ENABLE_REPORTING = false;

    /**
     * Subdirectory under the ant build directory where we store integration test md5 results
     */
    public static final String MD5_FILE_DB_SUBDIR = "integrationtests";

    public void setOutputFileLocation(File outputFileLocation) {
        this.outputFileLocation = outputFileLocation;
    }

    private static void ensureMd5DbDirectory() {
        // todo -- make path
        File dir = new File(MD5_FILE_DB_SUBDIR);
        if ( ! dir.exists() ) {
            System.out.printf("##### Creating MD5 db %s%n", MD5_FILE_DB_SUBDIR);
            if ( ! dir.mkdir() ) {
                throw new ReviewedStingException("Infrastructure failure: failed to create md5 directory " + MD5_FILE_DB_SUBDIR);
            }
        }
    }

    private static File getFileForMD5(final String md5) {
        final String basename = String.format("%s.integrationtest", md5);
        return new File(MD5_FILE_DB_SUBDIR + "/" + basename);
    }

    private static void updateMD5Db(final String md5, final File resultsFile) {
        // todo -- copy results file to DB dir if needed under filename for md5
        final File dbFile = getFileForMD5(md5);
        if ( ! dbFile.exists() ) {
            // the file isn't already in the db, copy it over
            System.out.printf("##### Updating MD5 file: %s%n", dbFile.getPath());
            try {
                FileUtils.copyFile(resultsFile, dbFile);
            } catch ( IOException e ) {
                throw new ReviewedStingException(e.getMessage());
            }
        } else {
            System.out.printf("##### MD5 file is up to date: %s%n", dbFile.getPath());

        }
    }

    private static String getMD5Path(final String md5, final String valueIfNotFound) {
        // todo -- look up the result in the directory and return the path if it exists
        final File dbFile = getFileForMD5(md5);
        return dbFile.exists() ? dbFile.getPath() : valueIfNotFound;
    }

    public String assertMatchingMD5(final String name, final File resultsFile, final String expectedMD5) {
        try {
            byte[] bytesOfMessage = getBytesFromFile(resultsFile);
            byte[] thedigest = MessageDigest.getInstance("MD5").digest(bytesOfMessage);
            BigInteger bigInt = new BigInteger(1, thedigest);
            String filemd5sum = bigInt.toString(16);
            while (filemd5sum.length() < 32) filemd5sum = "0" + filemd5sum; // pad to length 32

            //
            // copy md5 to integrationtests
            //
            updateMD5Db(filemd5sum, resultsFile);

            if (parameterize() || expectedMD5.equals("")) {
                System.out.println(String.format("PARAMETERIZATION[%s]: file %s has md5 = %s, stated expectation is %s, equal? = %b",
                                                 name, resultsFile, filemd5sum, expectedMD5, filemd5sum.equals(expectedMD5)));
            } else {
                System.out.println(String.format("Checking MD5 for %s [calculated=%s, expected=%s]", resultsFile, filemd5sum, expectedMD5));
                System.out.flush();

                if ( ! expectedMD5.equals(filemd5sum) ) {
                    // we are going to fail for real in assertEquals (so we are counted by the testing framework).
                    // prepare ourselves for the comparison
                    System.out.printf("##### Test %s is going fail #####%n", name);
                    String pathToExpectedMD5File = getMD5Path(expectedMD5, "[No DB file found]");
                    String pathToFileMD5File = getMD5Path(filemd5sum, "[No DB file found]");
                    System.out.printf("##### Path to expected   file (MD5=%s): %s%n", expectedMD5, pathToExpectedMD5File);
                    System.out.printf("##### Path to calculated file (MD5=%s): %s%n", filemd5sum, pathToFileMD5File);
                    System.out.printf("##### Diff command: diff %s %s%n", pathToExpectedMD5File, pathToFileMD5File);

                    // todo -- add support for simple inline display of the first N differences for text file
                }

                Assert.assertEquals(name + " Mismatching MD5s", expectedMD5, filemd5sum);
                System.out.println(String.format("  => %s PASSED", name));
            }

            return filemd5sum;
        } catch (Exception e) {
            throw new RuntimeException("Failed to read bytes from calls file: " + resultsFile, e);
        }
    }

    public List<String> assertMatchingMD5s(final String name, List<File> resultFiles, List<String> expectedMD5s) {
        List<String> md5s = new ArrayList<String>();
        for (int i = 0; i < resultFiles.size(); i++) {
            String md5 = assertMatchingMD5(name, resultFiles.get(i), expectedMD5s.get(i));
            md5s.add(i, md5);
        }

        return md5s;
    }


    public static byte[] getBytesFromFile(File file) throws IOException {
        InputStream is = new FileInputStream(file);

        // Get the size of the file
        long length = file.length();

        if (length > Integer.MAX_VALUE) {
            // File is too large
        }

        // Create the byte array to hold the data
        byte[] bytes = new byte[(int) length];

        // Read in the bytes
        int offset = 0;
        int numRead = 0;
        while (offset < bytes.length
                && (numRead = is.read(bytes, offset, bytes.length - offset)) >= 0) {
            offset += numRead;
        }

        // Ensure all the bytes have been read in
        if (offset < bytes.length) {
            throw new IOException("Could not completely read file " + file.getName());
        }

        // Close the input stream and return bytes
        is.close();
        return bytes;
    }

    public class WalkerTestSpec {
        String args = "";
        int nOutputFiles = -1;
        List<String> md5s = null;
        List<String> exts = null;
        Class expectedException = null;

        protected Map<String, File> auxillaryFiles = new HashMap<String, File>();

        public WalkerTestSpec(String args, int nOutputFiles, List<String> md5s) {
            this.args = args;
            this.nOutputFiles = nOutputFiles;
            this.md5s = md5s;
        }

        public WalkerTestSpec(String args, int nOutputFiles, List<String> exts, List<String> md5s) {
            this.args = args;
            this.nOutputFiles = nOutputFiles;
            this.md5s = md5s;
            this.exts = exts;
        }

        public WalkerTestSpec(String args, int nOutputFiles, Class expectedException) {
            this.args = args;
            this.nOutputFiles = nOutputFiles;
            this.expectedException = expectedException;
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

    protected Pair<List<File>, List<String>> executeTest(final String name, WalkerTestSpec spec) {
        ensureMd5DbDirectory(); // ensure the md5 directory exists

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
            return executeTest(name, null, tmpFiles, args, spec.getExpectedException());
        } else {
            List<String> md5s = new LinkedList<String>();
            md5s.addAll(spec.md5s);

            // check to see if they included any auxillary files, if so add them to the list
            for (String md5 : spec.auxillaryFiles.keySet()) {
                md5s.add(md5);
                tmpFiles.add(spec.auxillaryFiles.get(md5));
            }
            return executeTest(name, md5s, tmpFiles, args, null);
        }
    }

    public File createTempFile(String name, String extension) {
        try {
            File fl = File.createTempFile(name, extension);
            fl.deleteOnExit();
            return fl;
        } catch (IOException ex) {
            throw new ReviewedStingException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

    /**
     * execute the test, given the following:
     * @param name     the name of the test
     * @param md5s     the list of md5s
     * @param tmpFiles the temp file corresponding to the md5 list
     * @param args     the argument list
     * @return a pair of file and string lists
     */
    private Pair<List<File>, List<String>> executeTest(String name, List<String> md5s, List<File> tmpFiles, String args, Class expectedException) {
        CommandLineGATK instance = new CommandLineGATK();
        String[] command;

        // special case for ' and " so we can allow expressions
        if (args.indexOf('\'') != -1)
            command = escapeExpressions(args, "'");
        else if (args.indexOf('\"') != -1)
            command = escapeExpressions(args, "\"");
        else
            command = args.split(" ");

        if (outputFileLocation != null) {
            String[] cmd2 = Arrays.copyOf(command, command.length + 2);
            cmd2[command.length] = "-o";
            cmd2[command.length + 1] = this.outputFileLocation.getAbsolutePath();
            command = cmd2;
        }

        // add the logging level to each of the integration test commands
        String[] cmd2 = Arrays.copyOf(command, command.length + 4);
        cmd2[command.length] = "-l";
        cmd2[command.length+1] = "WARN";
        cmd2[command.length+2] = "-et";
        cmd2[command.length+3] = ENABLE_REPORTING ? "STANDARD" : "NO_ET";

        // run the executable
        boolean gotAnException = false;
        try {
            System.out.println(String.format("Executing test %s with GATK arguments: %s", name, Utils.join(" ",cmd2)));
            CommandLineExecutable.start(instance, cmd2);
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
            return null;
        } else {
            if ( CommandLineExecutable.result != 0) {
                throw new RuntimeException("Error running the GATK with arguments: " + args);
            }

            // clean up some memory
            instance = null;
            cmd2 = null;

            // we need to check MD5s
            return new Pair<List<File>, List<String>>(tmpFiles, assertMatchingMD5s(name, tmpFiles, md5s));
        }
    }

    private static String[] escapeExpressions(String args, String delimiter) {
        String[] command = {};
        String[] split = args.split(delimiter);
        for (int i = 0; i < split.length - 1; i += 2) {
            command = Utils.concatArrays(command, split[i].trim().split(" "));
            command = Utils.concatArrays(command, new String[]{split[i + 1]});
        }
        return Utils.concatArrays(command, split[split.length - 1].trim().split(" "));
    }

    @Test
    public void testWalkerUnitTest() {
        //System.out.println("WalkerTest is just a framework");
    }
}
