package org.broadinstitute.sting.gatk.datasources.sample;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Sep 9, 2010
 * Time: 8:21:00 AM
 */
public class SampleTest extends BaseTest {

    static Sample sampleA;
    static Sample sampleA1;
    static Sample sampleB;
    static Sample sampleC;

    @BeforeClass
    public static void init() {
        sampleA = new Sample("sampleA");
        sampleA.setProperty("uniqueProperty", "uniqueValue");
        sampleA1 = new Sample("sampleA");
        sampleA1.setProperty("uniqueProperty", "uniqueValue");
        sampleB = new Sample("sampleB");
        sampleC = new Sample("sampleC");
        sampleC.setProperty("population", "pop1");
        sampleC.setProperty("gender", Sample.Gender.MALE);
    }

    /**
     * Testing equality
     */
    @Test()
    public void equalsTest() {
        Assert.assertTrue(sampleA.equals(sampleA1));
        Assert.assertFalse(sampleA == sampleA1);
        Assert.assertFalse(sampleA.equals(sampleB));
    }

    /**
     * And hash
     */
    @Test()
    public void basicHashTest() {
        Assert.assertFalse(sampleA.hashCode() == sampleB.hashCode());
        Assert.assertTrue(sampleA.hashCode() == sampleA1.hashCode());
    }

    /**
     * Now test the special getter methods
     */
    @Test()
    public void specialGettersTest() {
        Assert.assertTrue(sampleC.getId().equals("sampleC"));
        Assert.assertTrue(sampleC.getPopulation().equals("pop1"));
        Assert.assertTrue(sampleC.isMale());
        Assert.assertFalse(sampleA.isMale());   // sample A doesn't have a gender, so this should be false
    }

}                     