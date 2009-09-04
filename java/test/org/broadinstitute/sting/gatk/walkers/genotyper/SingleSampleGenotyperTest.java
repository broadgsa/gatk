package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.genotype.Variant;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.junit.Test;
import org.broadinstitute.sting.gatk.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.executive.Accumulator;

import java.io.*;
import java.util.Arrays;
import java.util.EnumMap;
import java.security.MessageDigest;
import java.math.BigInteger;

import junit.framework.Assert;

public class SingleSampleGenotyperTest extends BaseTest {
    private final static double DELTA = 1e-8;

    private static final String oneMB = "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam";
    private static final String fiveMB = "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_15_mb.SLX.bam";
    private static String gatkargs(final String bam) {
        final String GATKArgs =
                "<GATK-argument-collection>\n" +
                        "   <sam-files class=\"java.util.ArrayList\">\n" +
                        "      <file>%s</file>\n" +
                        "   </sam-files>\n" +
                        "   <reference-file>/broad/1KG/reference/human_b36_both.fasta</reference-file>\n" +
                        "   <analysis-name>SingleSampleGenotyper</analysis-name>\n" +
                        "</GATK-argument-collection>";
        return String.format(GATKArgs, bam);
    }

    private static class ExpectedResult {
        long nCalls;
        String md5sum;
        public ExpectedResult(long calls, String md5) {
            this.nCalls = calls;
            this.md5sum = md5;
        }
    }

    public static EnumMap<BaseMismatchModel, ExpectedResult> expectedOutput = new EnumMap<BaseMismatchModel, ExpectedResult>(BaseMismatchModel.class);

    static {
        expectedOutput.put(BaseMismatchModel.ONE_STATE, new ExpectedResult(983L, "d5404668e76f206055f03d97162ea6d9"));
        expectedOutput.put(BaseMismatchModel.THREE_STATE, new ExpectedResult(1055L, "46fb7b66da3dac341e9c342f751d74cd"));
        expectedOutput.put(BaseMismatchModel.EMPIRICAL, new ExpectedResult(1070L, "ea0be2fd074a6c824a0670ad5b3e0aca"));
    }

    //private static double callingTolerance = 0.2; // I'm willing to tolerate +/- 10% variation in call rate from that expected

    public static long nExpectedCalls(BaseMismatchModel model) {
        return expectedOutput.get(model).nCalls;
    }

//    public static double nExpectedCallsTolerance(long nBasesCalled) {
//        return nExpectedCalls(nBasesCalled) * callingTolerance;
//    }

    public static void assertGoodNumberOfCalls(final String name, BaseMismatchModel model, SingleSampleGenotyper.CallResult calls) {
        Assert.assertEquals(name, nExpectedCalls(model), calls.nConfidentCalls);
    }   

    public static void assertMatchingMD5(final String name, BaseMismatchModel model, final File resultsFile ) {
        try {
            byte[] bytesOfMessage = getBytesFromFile(resultsFile);
            byte[] thedigest = MessageDigest.getInstance("MD5").digest(bytesOfMessage);
            BigInteger bigInt = new BigInteger(1, thedigest);
            String filemd5sum = bigInt.toString(16);
            Assert.assertEquals(name + "Mismatching MD5s", expectedOutput.get(model).md5sum, filemd5sum);
        } catch ( Exception e ) {
            throw new RuntimeException("Failed to read bytes from calls file: " + resultsFile);
        }
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

    private void oneMBCallTest(BaseMismatchModel model) {
        logger.warn("Executing oneMBCallTest for " + model);

        InputStream stream = new ByteArrayInputStream(gatkargs(oneMB).getBytes());
        GATKArgumentCollection mycoll = GATKArgumentCollection.unmarshal(stream);
        mycoll.intervals = Arrays.asList("1:10,000,000-11,000,000");

        SingleSampleGenotyper SSG = new SingleSampleGenotyper();
        SSG.VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;
        SSG.LOD_THRESHOLD = 5.0;

        try {
            SSG.VARIANTS_FILE = File.createTempFile("SSGTempTmpCalls.geli.calls", ".tmp" );
        } catch (IOException ex) {
            System.err.println("Cannot create temp file: " + ex.getMessage());
        }
        SSG.baseModel = model;

        GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        System.out.printf("model is %s, SSG says %s%n", model, SSG.baseModel );
        Accumulator obj = (Accumulator)engine.execute(mycoll, SSG);
        SingleSampleGenotyper.CallResult calls  = (SingleSampleGenotyper.CallResult)obj.finishTraversal();
        logger.warn(String.format("SSG(%s): made %d calls at %d bases", SSG.baseModel, calls.nConfidentCalls, calls.nCalledBases));
        assertMatchingMD5("testBasic", model, SSG.VARIANTS_FILE);
        assertGoodNumberOfCalls("testBasic", model, calls);
    }

    @Test public void oneStateOneMBTest() { oneMBCallTest(BaseMismatchModel.ONE_STATE); }
    @Test public void threeStateOneMBTest() { oneMBCallTest(BaseMismatchModel.THREE_STATE); }
    @Test public void empiricalStateOneMBTest() { oneMBCallTest(BaseMismatchModel.EMPIRICAL); }
}