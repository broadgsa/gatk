package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;

import org.testng.annotations.Test;

import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Jun 30, 2010
 * Time: 3:32:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class VCFHeaderUnitTest extends BaseTest {

    private VCFHeader createHeader(String headerStr) {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = (VCFHeader)codec.readHeader(new AsciiLineReader(new StringBufferInputStream(headerStr)));
        Assert.assertEquals(header.getMetaData().size(), VCF4headerStringCount);
        return header;
    }

    @Test
    public void testVCF4ToVCF4() {
        VCFHeader header = createHeader(VCF4headerStrings);
        checkMD5ofHeaderFile(header, "4648aa1169257e0a8a9d30131adb5f35");
    }

    @Test
    public void testVCF4ToVCF4_alternate() {
        VCFHeader header = createHeader(VCF4headerStrings_with_negativeOne);
        checkMD5ofHeaderFile(header, "ad8c4cf85e868b0261ab49ee2c613088");
    }

        /**
     * a little utility function for all tests to md5sum a file
     * Shameless taken from:
     *
     * http://www.javalobby.org/java/forums/t84420.html
     *
     * @param file the file
     * @return a string
     */
    private static String md5SumFile(File file) {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("MD5");
        } catch (NoSuchAlgorithmException e) {
            throw new ReviewedStingException("Unable to find MD5 digest");
        }
        InputStream is;
        try {
            is = new FileInputStream(file);
        } catch (FileNotFoundException e) {
            throw new ReviewedStingException("Unable to open file " + file);
        }
        byte[] buffer = new byte[8192];
        int read;
        try {
            while ((read = is.read(buffer)) > 0) {
                digest.update(buffer, 0, read);
            }
            byte[] md5sum = digest.digest();
            BigInteger bigInt = new BigInteger(1, md5sum);
            return bigInt.toString(16);

        }
        catch (IOException e) {
            throw new ReviewedStingException("Unable to process file for MD5", e);
        }
        finally {
            try {
                is.close();
            }
            catch (IOException e) {
                throw new ReviewedStingException("Unable to close input stream for MD5 calculation", e);
            }
        }
    }

    private void checkMD5ofHeaderFile(VCFHeader header, String md5sum) {
        File myTempFile = null;
        PrintWriter pw = null;
        try {
            myTempFile = File.createTempFile("VCFHeader","vcf");
            myTempFile.deleteOnExit();
            pw = new PrintWriter(myTempFile);
        } catch (IOException e) {
            Assert.fail("Unable to make a temp file!");
        }
        for (VCFHeaderLine line : header.getMetaData())
            pw.println(line);
        pw.close();
        Assert.assertTrue(md5sum.equals(md5SumFile(myTempFile)));
    }

    public static int VCF4headerStringCount = 16;

    public static String VCF4headerStrings =
                "##fileformat=VCFv4.0\n"+
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##INFO=<ID=CA, Number=1, Type=String, Description=\"Pilot 1 callability mask\">\n"+
                "##INFO=<ID=HP, Number=1, Type=Integer, Description=\"Reference homopolymer tract length\">\n"+
                "##INFO=<ID=NS, Number=1, Type=Integer, Description=\"Number of samples with data\">\n"+
                "##INFO=<ID=DB, Number=0, Type=Flag, Description=\"dbSNP membership build 129 - type match and indel sequence length match within 25 bp\">\n"+
                "##INFO=<ID=NR, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on reverse strand\">\n"+
                "##INFO=<ID=NF, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on forward strand\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";


    public static String VCF4headerStrings_with_negativeOne =
                "##fileformat=VCFv4.0\n"+
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=YY, Number=., Type=Integer, Description=\"Some weird value that has lots of parameters\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##INFO=<ID=CA, Number=1, Type=String, Description=\"Pilot 1 callability mask\">\n"+
                "##INFO=<ID=HP, Number=1, Type=Integer, Description=\"Reference homopolymer tract length\">\n"+
                "##INFO=<ID=NS, Number=1, Type=Integer, Description=\"Number of samples with data\">\n"+
                "##INFO=<ID=DB, Number=0, Type=Flag, Description=\"dbSNP membership build 129 - type match and indel sequence length match within 25 bp\">\n"+
                "##INFO=<ID=NR, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on reverse strand\">\n"+
                "##INFO=<ID=NF, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on forward strand\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=TT, Number=., Type=Integer, Description=\"Lots of TTs\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

}
