package org.broadinstitute.sting.gatk.refdata;
/*
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.oneoffprojects.variantcontext.Allele;
import org.broadinstitute.sting.oneoffprojects.variantcontext.Genotype;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.exceptions.StingException;
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
 * User: chartl
 * Date: Jan 22, 2010
 * Time: 11:27:33 PM
 * To change this template use File | Settings | File Templates.
 *
public class PlinkRodUnitTest extends BaseTest {
    // todo :: get the isIndel() isInsertion() and isDeletion() tests working again -- this may require new
    // todo :: methods in the objects themselves
    private static IndexedFastaSequenceFile seq;

    @BeforeClass
    public static void beforeTests() {
        try {
            seq = new IndexedFastaSequenceFile(new File(b36KGReference));
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
        PlinkRod rod = new PlinkRod("test");
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

        String alleleStr1 = new String(snp1_individual3_alleles.get(0).getBases())+new String(snp1_individual3_alleles.get(1).getBases());
        String alleleStr2 = new String(snp3_individual2_alleles.get(0).getBases())+new String(snp3_individual2_alleles.get(1).getBases());

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
        PlinkRod rod = new PlinkRod("test");
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
        // Assert.assertTrue("That the second genotype of second variant is not a point mutation", ! genotypesInRod.get(1).get(1). );
        // Assert.assertTrue("That the second genotype of fourth variant is not a point mutation", ! genotypesInRod.get(3).get(1).isPointGenotype() );
        Assert.assertTrue("That the second genotype of fourth variant is homozygous", genotypesInRod.get(3).get(1).isHom());
        Assert.assertTrue("That the fourth genotype of fourth variant is heterozygous",genotypesInRod.get(3).get(3).isHet());
        Assert.assertEquals("That the reference deletion genotype has the correct string", "ATTTAT",genotypesInRod.get(3).get(2).getAlleles().get(0).getBases());
        Assert.assertEquals("That the insertion bases are correct","CTC",genotypesInRod.get(1).get(2).getAlleles().get(0).getBases());
        Assert.assertEquals("That the snp bases are correct","GC",new String(genotypesInRod.get(2).get(2).getAlleles().get(0).getBases())+new String(genotypesInRod.get(2).get(2).getAlleles().get(1).getBases()));
    }

    @Test
    public void testBinaryPedFileNoIndels() {
        PlinkRod rod = new PlinkRod("binaryTest1");
        try {
            rod.initialize(new File("/humgen/gsa-hpprojects/GATK/data/Validation_Data/test/plink_rod_test/binary_noindel_test.bed"));
        } catch (FileNotFoundException e) {
            throw new StingException("Test file for testBinaryPedFileNoIndels() could not be found",e);
        }

        // iterate through the ROD and get stuff
        ArrayList<GenomeLoc> lociInRod = new ArrayList<GenomeLoc>();
        ArrayList<ArrayList<Genotype>> genotypesInRod = new ArrayList<ArrayList<Genotype>>();
        ArrayList<ArrayList<String>> samplesInRod = new ArrayList<ArrayList<String>>();

        do {
            lociInRod.add(rod.getLocation());
            genotypesInRod.add(rod.getGenotypes());
            samplesInRod.add(rod.getVariantSampleNames());
        } while ( rod.parseLine(null,null) );

        List<String> expecLoc = Arrays.asList("1:123456","1:14327877","2:22074511","3:134787","3:178678","4:829645","4:5234132","12:1268713");

        for ( int i = 0; i < expecLoc.size(); i ++ ) {
            Assert.assertEquals("That locus "+(i+1)+" in the rod is correct", expecLoc.get(i), lociInRod.get(i).toString());
        }

        List<String> expecAlleles = Arrays.asList("AA","AA","AA","GG","GG","GG","AA","TA","TT","CC","CC","GC","TC","CC","TT",
                                                  "GG","GG","AG","TT","CC","CT","TG","GG","GG");
        List<Boolean> expecHet = Arrays.asList(false,false,false,false,false,false,false,true,false,false,false,true,true,false,
                    false,false,false,true,false,false,true,true,false,false);
        List<String> expecName = Arrays.asList("NA12878","NA12890","NA07000","NA12878","NA12890","NA07000","NA12878","NA12890","NA07000",
                    "NA12878","NA12890","NA07000","NA12878","NA12890","NA07000","NA12878","NA12890","NA07000","NA12878","NA12890","NA07000",
                "NA12878","NA12890","NA07000");
        int snpNo = 1;
        int indiv = 1;
        int alleleOffset = 0;
        for ( ArrayList<Genotype> snp : genotypesInRod ) {
            for ( Genotype gen : snp ) {
                String alStr = new String(gen.getAlleles().get(0).getBases())+new String(gen.getAlleles().get(1).getBases());
                Assert.assertEquals("That the allele of person "+indiv+" for snp "+snpNo+" is correct "+
                                    "(allele offset "+alleleOffset+")", expecAlleles.get(alleleOffset),alStr);
                Assert.assertEquals("That the genotype of person "+indiv+" for snp "+snpNo+" is properly set", expecHet.get(alleleOffset),gen.isHet());
                Assert.assertEquals("That the name of person "+indiv+" for snp "+snpNo+" is correct", expecName.get(alleleOffset),samplesInRod.get(snpNo-1).get(indiv-1));
                indiv++;
                alleleOffset++;
            }
            indiv = 1;
            snpNo++;
        }
    }

    @Test
    public void testIndelBinary() {
        PlinkRod rod = new PlinkRod("binaryTest2");
        try {
            rod.initialize(new File("/humgen/gsa-hpprojects/GATK/data/Validation_Data/test/plink_rod_test/binary_indel_test.bed"));
        } catch (FileNotFoundException e) {
            throw new StingException("Test file for testBinaryPedFileNoIndels() could not be found",e);
        }

        ArrayList<ArrayList<Genotype>> genotypesInRod = new ArrayList<ArrayList<Genotype>>();
        do {
            genotypesInRod.add(rod.getGenotypes());
        } while ( rod.parseLine(null,null) );

        List<String> expecAlleles = Arrays.asList("ACCA","","ACCAACCA","GGGG","GG","","AA","TA","00","","CCTCCT","CCT",
                    "TC","CC","TT","GG","GG","AG","","CTTGCTTG","CTTG","TG","GG","GG");
        List<Boolean> expecDeletion = Arrays.asList(false,false,false,false,false,false,false,false,false,true,false,true,
                    false,false,false,false,false,false,true,false,true,false,false,false);
        List<Boolean> expecInsertion = Arrays.asList(true,false,true,true,true,false,false,false,false,false,false,false,
                    false,false,false,false,false,false,false,false,false,false,false,false);

        int al = 0;
        for ( ArrayList<Genotype> snp : genotypesInRod ) {
            for ( Genotype gen : snp ) {
                String alStr = new String(gen.getAlleles().get(0).getBases())+new String(gen.getAlleles().get(1).getBases());
                Allele firstAl = gen.getAlleles().get(0);
                Allele secondAl = gen.getAlleles().get(1);
                // boolean isInsertion = ( firstAl.getType() == Allele.AlleleType.INSERTION || secondAl.getType() == Allele.AlleleType.INSERTION );
                // boolean isDeletion =  ( firstAl.getType() == Allele.AlleleType.DELETION || secondAl.getType() == Allele.AlleleType.DELETION );
                Assert.assertEquals("That the indel file has the correct alleles for genotype "+al,expecAlleles.get(al), alStr);
                // Assert.assertEquals("That the insertion property of genotype "+al+" is correct",expecInsertion.get(al),isInsertion);
                // Assert.assertEquals("That the deletion property of genotype "+al+" is correct", expecDeletion.get(al), isDeletion);
                al++;
            }
        }
    }
}*/
