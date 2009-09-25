package org.broadinstitute.sting;

import junit.framework.Assert;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.Utils;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;

public class WalkerTest extends BaseTest {
    public String assertMatchingMD5(final String name, final File resultsFile, final String expectedMD5 ) {
        try {
            byte[] bytesOfMessage = getBytesFromFile(resultsFile);
            byte[] thedigest = MessageDigest.getInstance("MD5").digest(bytesOfMessage);
            BigInteger bigInt = new BigInteger(1, thedigest);
            String filemd5sum = bigInt.toString(16);
            while (filemd5sum.length() < 32) filemd5sum = "0" + filemd5sum; // pad to length 32
            if ( parameterize() || expectedMD5.equals("") ) {
                logger.warn(String.format("PARAMETERIZATION[%s]: file %s has md5 = %s, stated expectation is %s, equal? = %b",
                        name, resultsFile, filemd5sum, expectedMD5, filemd5sum.equals(expectedMD5)));
            } else {
                logger.warn(String.format("Checking MD5 for %s [calculated=%s, expected=%s]", resultsFile, filemd5sum, expectedMD5));
                Assert.assertEquals(name + " Mismatching MD5s", expectedMD5, filemd5sum);
                logger.warn(String.format("  => %s PASSED", name));
            }

            return filemd5sum;
        } catch ( Exception e ) {
            throw new RuntimeException("Failed to read bytes from calls file: " + resultsFile);
        }
    }

    public List<String> assertMatchingMD5s(final String name, List<File> resultFiles, List<String> expectedMD5s ) {
        List<String> md5s = new ArrayList<String>();
        for ( int i = 0; i < resultFiles.size(); i++ ) {
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
        byte[] bytes = new byte[(int)length];

        // Read in the bytes
        int offset = 0;
        int numRead = 0;
        while (offset < bytes.length
               && (numRead=is.read(bytes, offset, bytes.length-offset)) >= 0) {
            offset += numRead;
        }

        // Ensure all the bytes have been read in
        if (offset < bytes.length) {
            throw new IOException("Could not completely read file "+file.getName());
        }

        // Close the input stream and return bytes
        is.close();
        return bytes;
    }

    public class WalkerTestSpec {
        String args = "";
        int nOutputFiles = -1;
        List<String> md5s = null;

        public WalkerTestSpec(String args, int nOutputFiles, List<String> md5s) {
            this.args = args;
            this.nOutputFiles = nOutputFiles;
            this.md5s = md5s;
        }
    }

    protected boolean parameterize() {
        return false;
    }

    protected Pair<List<File>, List<String>> executeTest(final String name, WalkerTestSpec spec) {
        List<File> tmpFiles = new ArrayList<File>();
        for ( int i = 0; i < spec.nOutputFiles; i++ ) {
            try {
                File fl = File.createTempFile(String.format("walktest.tmp_param.%d", i), ".tmp" );
                fl.deleteOnExit();
                tmpFiles.add( fl );
            } catch (IOException ex) {
                System.err.println("Cannot create temp file: " + ex.getMessage());
            }
        }

        final String args = String.format(spec.args, tmpFiles.toArray());
        logger.warn(Utils.dupString('-', 80));
        logger.warn(String.format("Executing test %s with GATK arguments: %s", name, args));

        CommandLineGATK instance = new CommandLineGATK();
        CommandLineExecutable.start(instance, args.split(" "));

        if ( CommandLineExecutable.result != 0 ) {
            throw new RuntimeException("Error running the GATK with arguments: " + args);
        }

        return new Pair<List<File>, List<String>>(tmpFiles, assertMatchingMD5s(name, tmpFiles, spec.md5s));
    }

    @Test
    public void testWalkerTest() {
        //logger.warn("WalkerTest is just a framework");
    }
}