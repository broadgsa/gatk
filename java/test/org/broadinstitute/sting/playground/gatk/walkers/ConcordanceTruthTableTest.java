package org.broadinstitute.sting.playground.gatk.walkers;

import org.junit.Test;
import org.junit.Assert;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.BasicGenotype;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.playground.gatk.walkers.varianteval.ConcordanceTruthTable;
import org.broadinstitute.sting.BaseTest;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 19, 2009
 * Time: 9:40:51 AM
 * To change this template use File | Settings | File Templates.
 */
public class ConcordanceTruthTableTest extends BaseTest {

    @Test
    public void testAlleleFrequencyCalculation() {

        ConcordanceTruthTable ctt = new ConcordanceTruthTable(3);
        // this will test the count of non-reference alleles at a T/G polymorphic site
        Genotype ref1 = new BasicGenotype(null,"GG",'G',30);
        Genotype ref2 = new BasicGenotype(null,"GG",'G',30);
        Genotype ref3 = new BasicGenotype(null,"GG",'G',30);
        Genotype het1 = new BasicGenotype(null,"GT",'G',32);
        Genotype het2 = new BasicGenotype(null,"GT",'G',28);
        Genotype hom1 = new BasicGenotype(null,"TT",'G',40);
        Genotype hom2 = new BasicGenotype(null,"TT",'G',27);

        List<Pair<Genotype,Genotype>> oneHom = new ArrayList<Pair<Genotype,Genotype>>(4);
        oneHom.add(new Pair<Genotype,Genotype>(ref1,null));
        oneHom.add(new Pair<Genotype,Genotype>(null,null));
        oneHom.add(new Pair<Genotype,Genotype>(ref2,null));
        oneHom.add(new Pair<Genotype,Genotype>(ref3,null));
        oneHom.add(new Pair<Genotype,Genotype>(hom2,null));

        List<Pair<Genotype,Genotype>> oneHet = new ArrayList<Pair<Genotype,Genotype>>(4);
        oneHet.add(new Pair<Genotype,Genotype>(ref1,null));
        oneHet.add(new Pair<Genotype,Genotype>(ref2,null));
        oneHet.add(new Pair<Genotype,Genotype>(ref3,null));
        oneHet.add(new Pair<Genotype,Genotype>(null,null));
        oneHet.add(new Pair<Genotype,Genotype>(het1,null));

        List<Pair<Genotype,Genotype>> twoHetOneHom = new ArrayList<Pair<Genotype,Genotype>>(5);
        twoHetOneHom.add(new Pair<Genotype,Genotype>(ref1,null));
        twoHetOneHom.add(new Pair<Genotype,Genotype>(ref2,null));
        twoHetOneHom.add(new Pair<Genotype,Genotype>(ref3,null));
        twoHetOneHom.add(new Pair<Genotype,Genotype>(het1,null));
        twoHetOneHom.add(new Pair<Genotype,Genotype>(het2,null));
        twoHetOneHom.add(new Pair<Genotype,Genotype>(hom1,null));

        List<Pair<Genotype,Genotype>> twoHetTwoHom = new ArrayList<Pair<Genotype,Genotype>>(7);
        twoHetTwoHom.add(new Pair<Genotype,Genotype>(ref1,null));
        twoHetTwoHom.add(new Pair<Genotype,Genotype>(ref2,null));
        twoHetTwoHom.add(new Pair<Genotype,Genotype>(null,null));
        twoHetTwoHom.add(new Pair<Genotype,Genotype>(ref3,null));
        twoHetTwoHom.add(new Pair<Genotype,Genotype>(het1,null));
        twoHetTwoHom.add(new Pair<Genotype,Genotype>(het2,null));
        twoHetTwoHom.add(new Pair<Genotype,Genotype>(hom1,null));
        twoHetTwoHom.add(new Pair<Genotype,Genotype>(hom2,null));

        List<Pair<Genotype,Genotype>> hetHomNoRef = new ArrayList<Pair<Genotype,Genotype>>(2);
        hetHomNoRef.add(new Pair<Genotype,Genotype>(het2,null));
        hetHomNoRef.add(new Pair<Genotype,Genotype>(hom2,null));

        List<Pair<Genotype,Genotype>> homNoRef = new ArrayList<Pair<Genotype,Genotype>>(1);
        homNoRef.add(new Pair<Genotype,Genotype>(hom1,null));

        Pair<Genotype,Pair<Integer,Integer>> countShouldBeOne = ctt.getPooledAlleleFrequency(oneHet,'G');
        Pair<Genotype,Pair<Integer,Integer>>  countShouldBeTwo = ctt.getPooledAlleleFrequency(oneHom,'G');
        Pair<Genotype,Pair<Integer,Integer>>  countShouldBeFour = ctt.getPooledAlleleFrequency(twoHetOneHom,'G');
        Pair<Genotype,Pair<Integer,Integer>>  countShouldBeSix = ctt.getPooledAlleleFrequency(twoHetTwoHom,'G');
        Pair<Genotype,Pair<Integer,Integer>>  countShouldBeThree = ctt.getPooledAlleleFrequency(hetHomNoRef,'G');
        Pair<Genotype,Pair<Integer,Integer>>  countShouldBeTwoHereToo = ctt.getPooledAlleleFrequency(homNoRef, 'G');

        int expecChips = 4+4+6+7+2+1;
        int numChips = countShouldBeOne.getSecond().getSecond() + countShouldBeTwo.getSecond().getSecond() +
                       countShouldBeFour.getSecond().getSecond() + countShouldBeSix.getSecond().getSecond() +
                       countShouldBeThree.getSecond().getSecond() + countShouldBeTwoHereToo.getSecond().getSecond();


        logger.warn("Testing single het");
        Assert.assertTrue(countShouldBeOne.getSecond().getFirst() == 1);
        logger.warn("Testing single hom");
        Assert.assertTrue(countShouldBeTwo.getSecond().getFirst() == 2);
        logger.warn("Testing two hets + hom");
        Assert.assertTrue(countShouldBeFour.getSecond().getFirst() == 4);
        logger.warn("Testing two hets + two homs");
        Assert.assertTrue(countShouldBeSix.getSecond().getFirst() == 6);
        logger.warn("Testing het + hom without ref");
        Assert.assertTrue(countShouldBeThree.getSecond().getFirst() == 3);
        logger.warn("Testing hom without ref");
        Assert.assertTrue(countShouldBeTwoHereToo.getSecond().getFirst() == 2);
        logger.warn("Testing chip count sum");
        Assert.assertTrue( expecChips == numChips);
    }
}
