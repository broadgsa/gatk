package org.broadinstitute.sting.gatk.samples;

import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Sep 9, 2010
 * Time: 8:21:00 AM
 */
public class SampleUnitTest extends BaseTest {
    SampleDataSource db;
    static Sample fam1A, fam1B, fam1C;
    static Sample s1, s2;
    static Sample trait1, trait2, trait3, trait4;

    @BeforeClass
    public void init() {
        db = new SampleDataSource();

        fam1A = new Sample("1A", db, "fam1", "1B", "1C", Gender.UNKNOWN);
        fam1B = new Sample("1B", db, "fam1", null, null, Gender.MALE);
        fam1C = new Sample("1C", db, "fam1", null, null, Gender.FEMALE);

        s1 = new Sample("s1", db);
        s2 = new Sample("s2", db);

        trait1 = new Sample("t1", db, Affection.AFFECTED, Sample.UNSET_QT);
        trait2 = new Sample("t2", db, Affection.UNAFFECTED, Sample.UNSET_QT);
        trait3 = new Sample("t3", db, Affection.UNKNOWN, Sample.UNSET_QT);
        trait4 = new Sample("t4", db, Affection.QUANTITATIVE, 1.0);
    }

    /**
     * Now test the special getter methods
     */
    @Test()
    public void specialGettersTest() {
        // todo -- test for sample with extra properties, like population
//        Assert.assertTrue(sampleC.getID().equals("sampleC"));
//        Assert.assertTrue(sampleC.getPopulation().equals("pop1"));
    }

    @Test()
    public void testGenders() {
        Assert.assertTrue(fam1A.getGender() == Gender.UNKNOWN);
        Assert.assertTrue(fam1B.getGender() == Gender.MALE);
        Assert.assertTrue(fam1C.getGender() == Gender.FEMALE);
    }
}
