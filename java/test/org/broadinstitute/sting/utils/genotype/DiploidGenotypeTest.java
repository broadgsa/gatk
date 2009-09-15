package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Test;


/**
 * 
 * @author aaron 
 * 
 * Class DiploidGenotypeTest
 *
 * Testing the basic functionality of the diploid genotype class
 */
public class DiploidGenotypeTest extends BaseTest {

    @Test
    public void testCreateDiploidFromString() {
        String genotype = "AA";
        DiploidGenotype g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "AC";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "AG";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "AT";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "CC";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "CG";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "CT";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "GG";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "GT";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));

        genotype = "TT";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(genotype.equals(g.toString()));
    }

    @Test
    public void testIsHet() {
        String genotype = "AG";
        DiploidGenotype g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(g.isHet());

        genotype = "AA";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(!g.isHet());
    }

     @Test
    public void testIsHom() {
        String genotype = "AA";
        DiploidGenotype g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(g.isHom());

        genotype = "AG";
        g = DiploidGenotype.valueOf(genotype);
        Assert.assertTrue(!g.isHom());
    }

    @Test
      public void testCreateGenotype() {
        char ref = 'A';
        DiploidGenotype g = DiploidGenotype.createGenotype(ref);
        Assert.assertTrue("AA".equals(g.toString()));

        ref = 'a';
        g = DiploidGenotype.createGenotype(ref);
        Assert.assertTrue("AA".equals(g.toString()));

        ref = 't';
        g = DiploidGenotype.createGenotype(ref);
        Assert.assertTrue("TT".equals(g.toString()));

        ref = 'T';
        g = DiploidGenotype.createGenotype(ref);
        Assert.assertTrue("TT".equals(g.toString()));

    }


}
