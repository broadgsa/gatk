package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Assert;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: Jan 22, 2010
 * Time: 11:27:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class PlinkRodTest extends BaseTest {
    private static IndexedFastaSequenceFile seq;

    @BeforeClass
    public static void beforeTests() {
        try {
            seq = new IndexedFastaSequenceFile(new File(oneKGLocation + "reference/human_b36_both.fasta"));
        } catch (FileNotFoundException e) {
            throw new StingException("unable to load the sequence dictionary");
        }
        GenomeLocParser.setupRefContigOrdering(seq);

    }

    public BufferedReader openFile(String filename) {
        try {
            return new BufferedReader(new FileReader(filename));
        } catch (FileNotFoundException e) {
            throw new StingException("Couldn't open file " + filename);
        }

    }

    @Test
    public void testStandardPedFile() {
        PlinkRodWithGenomeLoc rod = new PlinkRodWithGenomeLoc("test");
        try {
            rod.initialize( new File("/humgen/gsa-hpprojects/GATK/data/Validation_Data/test/plink_rod_test/standard_plink_test.ped") );
        } catch ( FileNotFoundException e ) {
            throw new StingException("test file for testStandardPedFile() does not exist",e);
        }

        // test that the sample names are correct

        List<String> rodNames = rod.getVariantSampleNames();
        List<String> expectedNames = Arrays.asList("NA12887","NAMELY","COWBA");

        Assert.assertEquals("That there are as many samples in the rod as in the expected list",expectedNames.size(),rodNames.size());

        boolean namesCorrect = true;
        for ( int i = 0; i < expectedNames.size(); i++ ) {
            namesCorrect = namesCorrect && ( rodNames.get(i).equals(expectedNames.get(i)) );
        }

        Assert.assertTrue("That the names are correct and in the proper order",namesCorrect);

        // test that rod can be iterated over

        ArrayList<ArrayList<Genotype>> genotypesInRod = new ArrayList<ArrayList<Genotype>>();
        ArrayList<ArrayList<String>> sampleNamesInRod = new ArrayList<ArrayList<String>>();
        ArrayList<GenomeLoc> lociInRod = new ArrayList<GenomeLoc>();
        do {
            genotypesInRod.add(rod.getGenotypes());
            sampleNamesInRod.add(rod.getVariantSampleNames());
            lociInRod.add(rod.getLocation());
        } while ( rod.parseLine(null,null) );

        Assert.assertEquals("That there are three SNPs in the ROD",3,genotypesInRod.size());

        ArrayList<Genotype> snp1 = genotypesInRod.get(0);
        ArrayList<Genotype> snp3 = genotypesInRod.get(2);

        Assert.assertEquals("That there are three Genotypes in SNP 1",3,snp1.size());
        Assert.assertEquals("That there are three samples in SNP 2", 3, sampleNamesInRod.get(1).size());
        Assert.assertEquals("That there are three Genotypes in SNP 3",3,snp3.size());


        List<Allele> snp1_individual3_alleles = snp1.get(2).getAlleles();
        List<Allele> snp3_individual2_alleles = snp3.get(1).getAlleles();

        String alleleStr1 = snp1_individual3_alleles.get(0).getBases()+snp1_individual3_alleles.get(1).getBases();
        String alleleStr2 = snp3_individual2_alleles.get(0).getBases()+snp3_individual2_alleles.get(1).getBases();

        Assert.assertEquals("That the third genotype of snp 1 is correctly no-call","00",alleleStr1);
        Assert.assertEquals("That the second genotype of snp 3 is correctly G G","GG",alleleStr2);

        boolean name2isSame = true;

        for ( ArrayList<String> names : sampleNamesInRod ) {
            name2isSame = name2isSame && names.get(1).equals("NAMELY");
        }

        Assert.assertTrue("That the second name of all the genotypes is the same and is correct",name2isSame);

        // test that the loci are correctly parsed and in order

        List<String> expectedLoci = Arrays.asList("1:123456","2:13274","3:11111");
        boolean lociCorrect = true;
        for ( int i = 0; i < 3; i ++ ) {
            lociCorrect = lociCorrect && lociInRod.get(i).toString().equals(expectedLoci.get(i));
        }
    }

    @Test
    public void testStandardPedFileWithIndels() {
        PlinkRodWithGenomeLoc rod = new PlinkRodWithGenomeLoc("test");
        try {
            rod.initialize(new File("/humgen/gsa-hpprojects/GATK/data/Validation_Data/test/plink_rod_test/standard_plink_with_indels.ped") );
        } catch ( FileNotFoundException e) {
            throw new StingException("Test file for testStandardPedFileWithIndels() could not be found", e);
        }

        // Iterate through the rod

        List<ArrayList<Genotype>> genotypesInRod = new ArrayList<ArrayList<Genotype>>();
        ArrayList<ArrayList<String>> sampleNamesInRod = new ArrayList<ArrayList<String>>();
        ArrayList<GenomeLoc> lociInRod = new ArrayList<GenomeLoc>();
        ArrayList<Boolean> snpSites = new ArrayList<Boolean>();
        do {
            genotypesInRod.add(rod.getGenotypes());
            sampleNamesInRod.add(rod.getVariantSampleNames());
            lociInRod.add(rod.getLocation());
            snpSites.add(rod.variantIsSNP());
        } while ( rod.parseLine(null,null) );

        boolean snpOrder = true;
        List<Boolean> expectedOrder = Arrays.asList(true,false,true,false);
        for ( int i = 0; i < 4; i ++ ) {
            snpOrder = snpOrder && ( expectedOrder.get(i) == snpSites.get(i) );
        }

        Assert.assertTrue("That the variant type order is as expected", snpOrder);
        Assert.assertTrue("That the second genotype of second variant is not a point mutation", ! genotypesInRod.get(1).get(1).isPointGenotype() );
        Assert.assertTrue("That the second genotype of fourth variant is not a point mutation", ! genotypesInRod.get(3).get(1).isPointGenotype() );
        Assert.assertTrue("That the second genotype of fourth variant is homozygous", genotypesInRod.get(3).get(1).isHom());
        Assert.assertTrue("That the fourth genotype of fourth variant is heterozygous",genotypesInRod.get(3).get(3).isHet());
        Assert.assertEquals("That the reference deletion genotype has the correct string", "ATTTAT",genotypesInRod.get(3).get(2).getAlleles().get(0).getBases());
        Assert.assertEquals("That the insertion bases are correct","CTC",genotypesInRod.get(1).get(2).getAlleles().get(0).getBases());
        Assert.assertEquals("That the snp bases are correct","GC",genotypesInRod.get(2).get(2).getAlleles().get(0).getBases()+genotypesInRod.get(2).get(2).getAlleles().get(1).getBases());
    }

    @Test
    public void testBinaryPedFileNoIndels() {
        PlinkRodWithGenomeLoc rod = new PlinkRodWithGenomeLoc("binaryTest1");
        try {
            rod.initialize(new File("/humgen/gsa-hpprojects/GATK/data/Validation_Data/test/plink_rod_test/binary_noindel_test.bed"));
        } catch (FileNotFoundException e) {
            throw new StingException("Test file for testBinaryPedFileNoIndels() could not be found",e);
        }
    }
}
