// our package
package org.broadinstitute.sting.gatk.contexts.variantcontext;


// the imports for unit testing.

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ExpandingArrayList;
import org.broadinstitute.sting.utils.StingException;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

//    public Allele(byte[] bases, boolean isRef) {
//    public Allele(boolean isRef) {
//    public Allele(String bases, boolean isRef) {
//    public boolean isNullAllele()       { return length() == 0; }
//    public boolean isNonNullAllele()    { return ! isNullAllele(); }
//    public boolean isReference()        { return isRef; }
//    public boolean isNonReference()     { return ! isReference(); }
//    public byte[] getBases() { return bases; }
//    public boolean equals(Allele other) {
//    public int length() {

/**
 * Basic unit test for RecalData
 */
public class AlleleUnitTest extends BaseTest {
    Allele ARef, del, delRef, A, T, ATIns, ATCIns, NoCall;
    
    @Before
    public void before() {
        del = new Allele("-");
        delRef = new Allele("-", true);

        A = new Allele("A");
        ARef = new Allele("A", true);
        T = new Allele("T");

        ATIns = new Allele("AT");
        ATCIns = new Allele("ATC");

        NoCall = new Allele(".");
    }

    @Test
    public void testCreatingSNPAlleles() {
        logger.warn("testCreatingSNPAlleles");

        Assert.assertTrue(A.isNonReference());
        Assert.assertFalse(A.isReference());
        Assert.assertTrue(A.basesMatch("A"));
        Assert.assertEquals(A.length(), 1);
        Assert.assertTrue(A.isNonNull());
        Assert.assertFalse(A.isNull());

        Assert.assertTrue(ARef.isReference());
        Assert.assertFalse(ARef.isNonReference());
        Assert.assertTrue(ARef.basesMatch("A"));
        Assert.assertFalse(ARef.basesMatch("T"));

        Assert.assertTrue(T.isNonReference());
        Assert.assertFalse(T.isReference());
        Assert.assertTrue(T.basesMatch("T"));
        Assert.assertFalse(T.basesMatch("A"));
    }

    @Test
    public void testCreatingNoCallAlleles() {
        logger.warn("testCreatingNoCallAlleles");

        Assert.assertTrue(NoCall.isNonReference());
        Assert.assertFalse(NoCall.isReference());
        Assert.assertFalse(NoCall.basesMatch("."));
        Assert.assertEquals(NoCall.length(), 0);
        Assert.assertTrue(NoCall.isNonNull());
        Assert.assertFalse(NoCall.isNull());
    }


    @Test
    public void testCreatingIndelAlleles() {
        logger.warn("testCreatingIndelAlleles");

        Assert.assertEquals(ATIns.length(), 2);
        Assert.assertEquals(ATCIns.length(), 3);
        Assert.assertArrayEquals(ATIns.getBases(), "AT".getBytes());
        Assert.assertArrayEquals(ATCIns.getBases(), "ATC".getBytes());

        Assert.assertTrue(del.isNonReference());
        Assert.assertFalse(delRef.isNonReference());
        Assert.assertFalse(del.isReference());
        Assert.assertTrue(delRef.isReference());
        Assert.assertFalse(del.basesMatch("-"));
        Assert.assertTrue(del.basesMatch(""));
        Assert.assertEquals(del.length(), 0);
        Assert.assertFalse(del.isNonNull());
        Assert.assertTrue(del.isNull());
    }


    @Test
    public void testConstructors1() {
        logger.warn("testConstructors1");

        Allele a1 = new Allele("A");
        Allele a2 = new Allele("A".getBytes());
        Allele a3 = new Allele("a");
        Allele a4 = new Allele("A", true);

        Assert.assertTrue(a1.equals(a2));
        Assert.assertTrue(a1.equals(a3));
        Assert.assertFalse(a1.equals(a4));
    }

    @Test
    public void testDelConstructors() {
        logger.warn("testDelConstructors");

        Allele a1 = new Allele("-");
        Allele a2 = new Allele("-".getBytes());
        Allele a3 = new Allele("");
        Allele a4 = new Allele("", true);

        Assert.assertTrue(a1.equals(a2));
        Assert.assertTrue(a1.equals(a3));
        Assert.assertFalse(a1.equals(a4));
    }

    @Test
    public void testInsConstructors() {
        logger.warn("testInsConstructors");

        Allele a1 = new Allele("AC");
        Allele a2 = new Allele("AC".getBytes());
        Allele a3 = new Allele("Ac");
        Allele a4 = new Allele("AC", true);

        Assert.assertTrue(a1.equals(a2));
        Assert.assertTrue(a1.equals(a3));
        Assert.assertFalse(a1.equals(a4));
    }

    @Test
    public void testEquals() {
        logger.warn("testEquals");
        Assert.assertTrue(ARef.basesMatch(A));
        Assert.assertFalse(ARef.equals(A));
        Assert.assertFalse(ARef.equals(del));
        Assert.assertFalse(ARef.equals(ATIns));
        Assert.assertFalse(ARef.equals(ATCIns));

        Assert.assertTrue(T.basesMatch(T));
        Assert.assertFalse(T.basesMatch(A));
        Assert.assertFalse(T.equals(A));

        Assert.assertTrue(del.basesMatch(del));
        Assert.assertTrue(del.basesMatch(delRef));
        Assert.assertTrue(del.equals(del));
        Assert.assertFalse(del.equals(delRef));

        Assert.assertTrue(ATIns.equals(ATIns));
        Assert.assertFalse(ATIns.equals(ATCIns));
        Assert.assertTrue(ATIns.basesMatch("AT"));
        Assert.assertFalse(ATIns.basesMatch("A"));
        Assert.assertFalse(ATIns.basesMatch("ATC"));

        Assert.assertTrue(ATIns.basesMatch("at"));
        Assert.assertFalse(ATIns.basesMatch("atc"));
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgs1() {
        logger.warn("testBadConstructorArgs1");
        byte[] foo = null;
        new Allele(foo);
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgs2() {
        logger.warn("testBadConstructorArgs2");
        new Allele("x");
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgs3() {
        logger.warn("testBadConstructorArgs3");
        new Allele("--");
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgs4() {
        logger.warn("testBadConstructorArgs4");
        new Allele("-A");
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgs5() {
        logger.warn("testBadConstructorArgs5");
        new Allele("A A");
    }
}