package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Test;


/**
 * @author aaron
 *         <p/>
 *         Class VCFGenotypeEncodingUnitTest
 *         <p/>
 *         test the VCFGenotypeEncoding class
 */
public class VCFGenotypeEncodingUnitTest extends BaseTest {
    @Test
    public void testDecodingSingle() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("A");
        Assert.assertTrue("A".equals(enc.toString()));
        Assert.assertEquals(0, enc.getLength());
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.SINGLE_BASE, enc.getType());

        VCFGenotypeEncoding enc2 = new VCFGenotypeEncoding("C");
        Assert.assertTrue("C".equals(enc2.toString()));
        Assert.assertEquals(0, enc.getLength());
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.SINGLE_BASE, enc.getType());

        VCFGenotypeEncoding enc3 = new VCFGenotypeEncoding("G");
        Assert.assertTrue("G".equals(enc3.toString()));
        Assert.assertEquals(0, enc.getLength());
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.SINGLE_BASE, enc.getType());

        VCFGenotypeEncoding enc4 = new VCFGenotypeEncoding("T");
        Assert.assertTrue("T".equals(enc4.toString()));
        Assert.assertEquals(0, enc.getLength());
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.SINGLE_BASE, enc.getType());

        VCFGenotypeEncoding enc5 = new VCFGenotypeEncoding("a");
        Assert.assertTrue("A".equals(enc5.toString()));
        Assert.assertEquals(0, enc.getLength());
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.SINGLE_BASE, enc.getType());

        VCFGenotypeEncoding enc6 = new VCFGenotypeEncoding("c");
        Assert.assertTrue("C".equals(enc6.toString()));
        Assert.assertEquals(0, enc.getLength());
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.SINGLE_BASE, enc.getType());

        VCFGenotypeEncoding enc7 = new VCFGenotypeEncoding("g");
        Assert.assertTrue("G".equals(enc7.toString()));
        Assert.assertEquals(0, enc.getLength());
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.SINGLE_BASE, enc.getType());

        VCFGenotypeEncoding enc8 = new VCFGenotypeEncoding("t");
        Assert.assertTrue("T".equals(enc8.toString()));
        Assert.assertEquals(0, enc.getLength());
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.SINGLE_BASE, enc.getType());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testDecodingSingleBadBase() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("E");
    }

    @Test(expected = IllegalArgumentException.class)
    public void testDecodingSingleWrongBase() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("I");
    }

    @Test
    public void testValidIndel() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("IAGGC");
        Assert.assertEquals(4, enc.getLength());
        Assert.assertTrue(enc.getBases().equals("AGGC"));
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.INSERTION, enc.getType());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testBadIndel() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("IAGRC");
    }

    @Test
    public void testValidDel() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("D40");
        Assert.assertEquals(40, enc.getLength());
        Assert.assertTrue(enc.getBases().equals(""));
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.DELETION, enc.getType());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testBadDel() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("DAGCT");
    }

    @Test
    public void testValidNoCall() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding(".");
        Assert.assertEquals(0, enc.getLength());
        Assert.assertTrue(enc.getBases().equals("."));
        Assert.assertEquals(VCFGenotypeEncoding.TYPE.UNCALLED, enc.getType());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testBadNoCall() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("..");
    }

    @Test
    public void testEquals() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("A");
        VCFGenotypeEncoding enc2 = new VCFGenotypeEncoding("A");
        VCFGenotypeEncoding enc3 = new VCFGenotypeEncoding("C");
        Assert.assertTrue(enc.equals(enc2));
        Assert.assertTrue(!enc.equals(enc3));
        enc = new VCFGenotypeEncoding("D40");
        enc2 = new VCFGenotypeEncoding("D40");
        enc3 = new VCFGenotypeEncoding("D41");
        Assert.assertTrue(enc.equals(enc2));
        Assert.assertTrue(!enc.equals(enc3));
        enc = new VCFGenotypeEncoding("IAAC");
        enc2 = new VCFGenotypeEncoding("IAAC");
        enc3 = new VCFGenotypeEncoding("IACG");
        Assert.assertTrue(enc.equals(enc2));
        Assert.assertTrue(!enc.equals(enc3));
        enc = new VCFGenotypeEncoding(".");
        enc2 = new VCFGenotypeEncoding(".");
        Assert.assertTrue(enc.equals(enc2));
    }

    @Test
    public void testHashCode() {
        VCFGenotypeEncoding enc = new VCFGenotypeEncoding("A");
        VCFGenotypeEncoding enc2 = new VCFGenotypeEncoding("A");
        VCFGenotypeEncoding enc3 = new VCFGenotypeEncoding("C");
        Assert.assertTrue(enc.hashCode() == enc2.hashCode());
        Assert.assertTrue(enc.hashCode() != enc3.hashCode());
        enc = new VCFGenotypeEncoding("D40");
        enc2 = new VCFGenotypeEncoding("D40");
        enc3 = new VCFGenotypeEncoding("D41");
        Assert.assertTrue(enc.hashCode() == enc2.hashCode());
        Assert.assertTrue(enc.hashCode() != enc3.hashCode());
        enc = new VCFGenotypeEncoding("IAAC");
        enc2 = new VCFGenotypeEncoding("IAAC");
        enc3 = new VCFGenotypeEncoding("IACG");
        Assert.assertTrue(enc.hashCode() == enc2.hashCode());
        Assert.assertTrue(enc.hashCode() != enc3.hashCode());
        enc = new VCFGenotypeEncoding(".");
        enc2 = new VCFGenotypeEncoding(".");
        Assert.assertTrue(enc.hashCode() == enc2.hashCode());
    }
}
