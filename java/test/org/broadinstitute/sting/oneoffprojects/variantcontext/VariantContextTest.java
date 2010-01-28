// our package
package org.broadinstitute.sting.oneoffprojects.variantcontext;


// the imports for unit testing.

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.junit.BeforeClass;

import java.util.Arrays;
import java.util.List;
import java.io.FileNotFoundException;
import java.io.File;

import net.sf.picard.reference.ReferenceSequenceFile;

/**
 * Basic unit test for RecalData
 */
public class VariantContextTest extends BaseTest {
    private static ReferenceSequenceFile seq;

    @BeforeClass
    public static void init() throws FileNotFoundException {
        // sequence
        seq = new IndexedFastaSequenceFile(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq);
    }
    
    Allele A, Aref, T, Tref;
    Allele del, delRef, ATC, ATCref;

    // A [ref] / T at 10
    GenomeLoc snpLoc = GenomeLocParser.createGenomeLoc("chr1", 10, 11);

    // - / ATC [ref] from 20-23
    GenomeLoc delLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 22);

    // - [ref] / ATC from 20-20
    GenomeLoc insLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 20);

    // - / A / T / ATC [ref] from 20-23
    GenomeLoc mixedLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 22);

    @Before
    public void before() {
        del = new Allele("-");
        delRef = new Allele("-", true);

        A = new Allele("A");
        Aref = new Allele("A", true);
        T = new Allele("T");
        Tref = new Allele("T", true);

        ATC = new Allele("ATC");
        ATCref = new Allele("ATC", true);
    }

    @Test
    public void testCreatingSNPVariantContext() {
        logger.warn("testCreatingSNPVariantContext");

        List<Allele> alleles = Arrays.asList(Aref, T);
        VariantContext vc = new VariantContext(snpLoc, alleles);
        logger.warn("vc = " + vc);

        Assert.assertEquals(vc.getLocation(), snpLoc);
        Assert.assertEquals(vc.getType(), VariantContext.Type.SNP);
        Assert.assertTrue(vc.isSNP());
        Assert.assertFalse(vc.isIndel());
        Assert.assertFalse(vc.isInsertion());
        Assert.assertFalse(vc.isDeletion());
        Assert.assertFalse(vc.isMixed());
        Assert.assertTrue(vc.isBiallelic());
        Assert.assertEquals(vc.getNAlleles(), 2);

        Assert.assertEquals(vc.getReference(), Aref);
        Assert.assertEquals(vc.getAlleles().size(), 2);
        Assert.assertEquals(vc.getAlternateAlleles().size(), 1);
        Assert.assertEquals(vc.getAlternateAllele(0), T);

        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingRefVariantContext() {
        logger.warn("testCreatingRefVariantContext");

        List<Allele> alleles = Arrays.asList(Aref);
        VariantContext vc = new VariantContext(snpLoc, alleles);
        logger.warn("vc = " + vc);

        Assert.assertEquals(vc.getLocation(), snpLoc);
        Assert.assertEquals(vc.getType(), VariantContext.Type.NO_VARIATION);
        Assert.assertFalse(vc.isSNP());
        Assert.assertFalse(vc.isIndel());
        Assert.assertFalse(vc.isInsertion());
        Assert.assertFalse(vc.isDeletion());
        Assert.assertFalse(vc.isMixed());
        Assert.assertFalse(vc.isBiallelic());
        Assert.assertEquals(vc.getNAlleles(), 1);

        Assert.assertEquals(vc.getReference(), Aref);
        Assert.assertEquals(vc.getAlleles().size(), 1);
        Assert.assertEquals(vc.getAlternateAlleles().size(), 0);
        //Assert.assertEquals(vc.getAlternateAllele(0), T);

        Assert.assertFalse(vc.hasGenotypes());
        Assert.assertEquals(vc.getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingDeletionVariantContext() {
        logger.warn("testCreatingDeletionVariantContext");

        List<Allele> alleles = Arrays.asList(ATCref, del);
        VariantContext vc = new VariantContext(delLoc, alleles);
        logger.warn("vc = " + vc);

        Assert.assertEquals(vc.getLocation(), delLoc);
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);
        Assert.assertFalse(vc.isSNP());
        Assert.assertTrue(vc.isIndel());
        Assert.assertFalse(vc.isInsertion());
        Assert.assertTrue(vc.isDeletion());
        Assert.assertFalse(vc.isMixed());
        Assert.assertTrue(vc.isBiallelic());
        Assert.assertEquals(vc.getNAlleles(), 2);

        Assert.assertEquals(vc.getReference(), ATCref);
        Assert.assertEquals(vc.getAlleles().size(), 2);
        Assert.assertEquals(vc.getAlternateAlleles().size(), 1);
        Assert.assertEquals(vc.getAlternateAllele(0), del);

        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingInsertionVariantContext() {
        logger.warn("testCreatingInsertionVariantContext");

        List<Allele> alleles = Arrays.asList(delRef, ATC);
        VariantContext vc = new VariantContext(insLoc, alleles);
        logger.warn("vc = " + vc);

        Assert.assertEquals(vc.getLocation(), insLoc);
        Assert.assertEquals(vc.getType(), VariantContext.Type.INDEL);
        Assert.assertFalse(vc.isSNP());
        Assert.assertTrue(vc.isIndel());
        Assert.assertTrue(vc.isInsertion());
        Assert.assertFalse(vc.isDeletion());
        Assert.assertFalse(vc.isMixed());
        Assert.assertTrue(vc.isBiallelic());
        Assert.assertEquals(vc.getNAlleles(), 2);

        Assert.assertEquals(vc.getReference(), delRef);
        Assert.assertEquals(vc.getAlleles().size(), 2);
        Assert.assertEquals(vc.getAlternateAlleles().size(), 1);
        Assert.assertEquals(vc.getAlternateAllele(0), ATC);

        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.getSampleNames().size(), 0);
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgs1() {
        logger.warn("testBadConstructorArgs1");
        new VariantContext(insLoc, Arrays.asList(delRef, ATCref));
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgs2() {
        logger.warn("testBadConstructorArgs2");
        new VariantContext(insLoc, Arrays.asList(delRef, del));
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgs3() {
        logger.warn("testBadConstructorArgs3");
        new VariantContext(insLoc, Arrays.asList(del));
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgsDuplicateAlleles1() {
        logger.warn("testBadConstructorArgsDuplicateAlleles1");
        new VariantContext(insLoc, Arrays.asList(Aref, T, T));
    }

    @Test (expected = IllegalArgumentException.class)
    public void testBadConstructorArgsDuplicateAlleles2() {
        logger.warn("testBadConstructorArgsDuplicateAlleles2");
        new VariantContext(insLoc, Arrays.asList(Aref, A));
    }


    @Test
    public void testAccessingSimpleSNPGenotypes() {
        logger.warn("testAccessingSimpleSNPGenotypes");

        List<Allele> alleles = Arrays.asList(Aref, T);
        VariantContext vc = new VariantContext(snpLoc, alleles);
        logger.warn("vc = " + vc);

        Genotype g1 = new Genotype(Arrays.asList(Aref, Aref), "AA", 10);
        Genotype g2 = new Genotype(Arrays.asList(Aref, T), "AT", 10);
        Genotype g3 = new Genotype(Arrays.asList(T, T), "TT", 10);

        vc.addGenotypes(Arrays.asList(g1, g2, g3));

        Assert.assertTrue(vc.hasGenotypes());
        Assert.assertFalse(vc.isMonomorphic());
        Assert.assertTrue(vc.isPolymorphic());
        Assert.assertEquals(vc.getSampleNames().size(), 3);

        Assert.assertEquals(vc.getGenotypes().size(), 3);
        Assert.assertEquals(vc.getGenotypes().get("AA"), g1);
        Assert.assertEquals(vc.getGenotype("AA"), g1);
        Assert.assertEquals(vc.getGenotypes().get("AT"), g2);
        Assert.assertEquals(vc.getGenotype("AT"), g2);
        Assert.assertEquals(vc.getGenotypes().get("TT"), g3);
        Assert.assertEquals(vc.getGenotype("TT"), g3);

        Assert.assertTrue(vc.hasGenotype("AA"));
        Assert.assertTrue(vc.hasGenotype("AT"));
        Assert.assertTrue(vc.hasGenotype("TT"));
        Assert.assertFalse(vc.hasGenotype("foo"));
        Assert.assertFalse(vc.hasGenotype("TTT"));
        Assert.assertFalse(vc.hasGenotype("at"));
        Assert.assertFalse(vc.hasGenotype("tt"));

        Assert.assertEquals(vc.getChromosomeCount(), 6);
        Assert.assertEquals(vc.getChromosomeCount(Aref), 3);
        Assert.assertEquals(vc.getChromosomeCount(T), 3);
    }

    @Test
    public void testAccessingCompleteGenotypes() {
        logger.warn("testAccessingCompleteGenotypes");

        List<Allele> alleles = Arrays.asList(Aref, T, del);
        VariantContext vc = new VariantContext(delLoc, alleles);
        logger.warn("vc = " + vc);

        Genotype g1 = new Genotype(Arrays.asList(Aref, Aref), "AA", 10);
        Genotype g2 = new Genotype(Arrays.asList(Aref, T), "AT", 10);
        Genotype g3 = new Genotype(Arrays.asList(T, T), "TT", 10);
        Genotype g4 = new Genotype(Arrays.asList(T, del), "Td", 10);
        Genotype g5 = new Genotype(Arrays.asList(del, del), "dd", 10);
        Genotype g6 = new Genotype(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), "..", 10);

        vc.addGenotypes(Arrays.asList(g1, g2, g3, g4, g5, g6));

        Assert.assertTrue(vc.hasGenotypes());
        Assert.assertFalse(vc.isMonomorphic());
        Assert.assertTrue(vc.isPolymorphic());
        Assert.assertEquals(vc.getGenotypes().size(), 6);

        Assert.assertEquals(10, vc.getChromosomeCount());
        Assert.assertEquals(3, vc.getChromosomeCount(Aref));
        Assert.assertEquals(4, vc.getChromosomeCount(T));
        Assert.assertEquals(3, vc.getChromosomeCount(del));
        Assert.assertEquals(2, vc.getChromosomeCount(Allele.NO_CALL));
    }

    @Test
    public void testAccessingRefGenotypes() {
        logger.warn("testAccessingRefGenotypes");

        List<Allele> alleles1 = Arrays.asList(Aref, T);
        List<Allele> alleles2 = Arrays.asList(Aref);
        List<Allele> alleles3 = Arrays.asList(Aref, T, del);
        for ( List<Allele> alleles : Arrays.asList(alleles1, alleles2, alleles3)) {
            VariantContext vc = new VariantContext(snpLoc, alleles);
            logger.warn("vc = " + vc);

            Genotype g1 = new Genotype(Arrays.asList(Aref, Aref), "AA1", 10);
            Genotype g2 = new Genotype(Arrays.asList(Aref, Aref), "AA2", 10);
            Genotype g3 = new Genotype(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), "..", 10);

            vc.addGenotypes(Arrays.asList(g1, g2, g3));

            Assert.assertTrue(vc.hasGenotypes());
            Assert.assertTrue(vc.isMonomorphic());
            Assert.assertFalse(vc.isPolymorphic());
            Assert.assertEquals(vc.getGenotypes().size(), 3);

            Assert.assertEquals(4, vc.getChromosomeCount());
            Assert.assertEquals(4, vc.getChromosomeCount(Aref));
            Assert.assertEquals(0, vc.getChromosomeCount(T));
            Assert.assertEquals(2, vc.getChromosomeCount(Allele.NO_CALL));
        }
    }


    @Test
    public void testManipulatingAlleles() {
        logger.warn("testManipulatingAlleles");
        // todo -- add tests that call add/set/remove
    }

    @Test
    public void testManipulatingGenotypes() {
        logger.warn("testManipulatingGenotypes");
        // todo -- add tests that call add/set/remove
    }
}

// genotype functions
//    public boolean hasGenotypes() { return genotypes.size() > 0; }
//    public Map<String, Genotype> getGenotypes() { return genotypes; }
//    public Set<String> getSampleNames() {
//    public int getChromosomeCount() {
//    public int getChromosomeCount(Allele a) {
//    public boolean isMonomorphic() {
//    public boolean isPolymorphic() {
//    public Genotype getGenotype(String sample) {
//    public boolean hasGenotype(String sample) {
//    public void setGenotypes(Genotype genotype) {
//    public void setGenotypes(Collection<Genotype> genotypes) {
//    public void setGenotypes(Map<String, Genotype> genotypes) {
//    public void addGenotype(Genotype genotype) {
//    public void addGenotype(String sampleName, Genotype genotype) {
//    public void addGenotype(String sampleName, Genotype genotype, boolean allowOverwrites) {
//    public void removeGenotype(String sampleName) {
//    public void removeGenotype(Genotype genotype) {

// all functions
//    public Type getType() {
//    public boolean isSNP() { return getType() == Type.SNP; }
//    public boolean isVariant() { return getType() != Type.NO_VARIATION; }
//    public boolean isIndel() { return getType() == Type.INDEL; }
//    public boolean isMixed() { return getType() == Type.MIXED; }
//    public GenomeLoc getLocation() { return loc; }
//    public Allele getReference() {
//    public boolean isBiallelic() {
//    public boolean isMonomorphic() {
//    public boolean isPolymorphic() {
//    public int getNAlleles() {
//    public Set<Allele> getAlleles() { return alleles; }
//    public Set<Allele> getAlternateAlleles() {
//    public Allele getAlternateAllele(int i) {
//    public void setAlleles(Set<Allele> alleles) {
//    public void addAllele(Allele allele) {
//    public void addAllele(Allele allele, boolean allowDuplicates) {

//    public boolean validate() {
