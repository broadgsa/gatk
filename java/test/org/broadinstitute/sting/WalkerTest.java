package org.broadinstitute.sting;

import junit.framework.Assert;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.Utils;
import org.junit.Test;
import org.apache.log4j.Appender;
import org.apache.log4j.WriterAppender;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.util.*;

public class WalkerTest extends BaseTest {
    public String assertMatchingMD5(final String name, final File resultsFile, final String expectedMD5 ) {
        try {
            byte[] bytesOfMessage = getBytesFromFile(resultsFile);
            byte[] thedigest = MessageDigest.getInstance("MD5").digest(bytesOfMessage);
            BigInteger bigInt = new BigInteger(1, thedigest);
            String filemd5sum = bigInt.toString(16);
            while (filemd5sum.length() < 32) filemd5sum = "0" + filemd5sum; // pad to length 32
            if ( parameterize() || expectedMD5.equals("") ) {
                System.out.println(String.format("PARAMETERIZATION[%s]: file %s has md5 = %s, stated expectation is %s, equal? = %b",
                        name, resultsFile, filemd5sum, expectedMD5, filemd5sum.equals(expectedMD5)));
            } else {
                System.out.println(String.format("Checking MD5 for %s [calculated=%s, expected=%s]", resultsFile, filemd5sum, expectedMD5));
                System.out.flush();
                Assert.assertEquals(name + " Mismatching MD5s", expectedMD5, filemd5sum);
                System.out.println(String.format("  => %s PASSED", name));
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
        List<String> exts = null;

        Map<String,File> auxillaryFiles = new HashMap<String,File>();

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

        public void addAuxFile(String expectededMD5sum, File outputfile) {
            auxillaryFiles.put(expectededMD5sum,outputfile);
        }
    }

    protected boolean parameterize() {
        return false;
    }

    protected Pair<List<File>, List<String>> executeTest(final String name, WalkerTestSpec spec) {
        List<File> tmpFiles = new ArrayList<File>();
        for ( int i = 0; i < spec.nOutputFiles; i++ ) {
            try {
                String ext = spec.exts == null ? ".tmp" : "." + spec.exts.get(i);
                File fl = File.createTempFile(String.format("walktest.tmp_param.%d", i), ext );
                fl.deleteOnExit();
                tmpFiles.add( fl );
            } catch (IOException ex) {
                System.err.println("Cannot create temp file: " + ex.getMessage());
            }
        }

        final String args = String.format(spec.args, tmpFiles.toArray());
        System.out.println(Utils.dupString('-', 80));
        System.out.println(String.format("Executing test %s with GATK arguments: %s", name, args));

        List<String> md5s = new LinkedList<String>();
        md5s.addAll(spec.md5s);

        // check to see if they included any auxillary files, if so add them to the list
        for (String md5 : spec.auxillaryFiles.keySet()) {
            md5s.add(md5);
            tmpFiles.add(spec.auxillaryFiles.get(md5));
        }

        return executeTest(name, md5s, tmpFiles, args);
    }


    /**
     * execute the test, given the following:
     * @param name the name of the test
     * @param md5s the list of md5s
     * @param tmpFiles the temp file corresponding to the md5 list
     * @param args the argument list
     * @return a pair of file and string lists
     */
    private Pair<List<File>, List<String>> executeTest(String name, List<String> md5s, List<File> tmpFiles, String args) {
        CommandLineGATK instance = new CommandLineGATK();
        String[] command;

        // special case for ' and " so we can allow expressions
        if ( args.indexOf('\'') != -1 )
            command = escapeExpressions(args, "'");
        else if ( args.indexOf('\"') != -1 )
            command = escapeExpressions(args, "\"");
        else
            command = args.split(" ");

        CommandLineExecutable.start(instance, command);

        if ( CommandLineExecutable.result != 0 ) {
            throw new RuntimeException("Error running the GATK with arguments: " + args);
        }

        return new Pair<List<File>, List<String>>(tmpFiles, assertMatchingMD5s(name, tmpFiles, md5s));
    }

    private static String[] escapeExpressions(String args, String delimiter) {
        String[] command = {};
        String[] split = args.split(delimiter);
        for (int i = 0; i < split.length - 1; i+=2) {
            command = Utils.concatArrays(command, split[i].trim().split(" "));
            command = Utils.concatArrays(command, new String[] { split[i+1] } );
        }
        return Utils.concatArrays(command, split[split.length-1].trim().split(" "));
    }

    @Test
    public void testWalkerTest() {
        //System.out.println("WalkerTest is just a framework");
    }
}