// our package
package org.broadinstitute.sting.gatk.walkers.recalibration;

import java.util.*;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;

// the imports for unit testing.
import org.junit.Assert;
import org.junit.Test;
import org.junit.Before;
import org.junit.BeforeClass;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.playground.gatk.walkers.indels.CleanedReadInjector;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMFileReader;
import org.broadinstitute.sting.utils.sam.ArtificialSAMFileWriter;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.HashSet;
import java.io.FileNotFoundException;
import java.io.File;

/**
 * Basic unit test for RecalData
 */
public class CovariateCounterTest extends BaseTest {
    String readGroup1 = "rg1";
    String readGroup2 = "rg2";
    Set<String> readGroups = new HashSet<String>();

    SAMFileHeader header;

    SAMRecord read1, read2, read3;

    byte bases1[] = {'a', 't', 'c', 'g', 'a'};
    byte quals1[] = {1, 2, 3, 4, 5};
    byte quals3[] = {1, 2, 5, 5, 5};
    byte bases2[] = {'t', 'c', 'g', 'a', 't'};
    byte quals2[] = {2, 2, 4, 5, 2};

    /*
    public CovariateCounter( Set<String> readGroups, boolean collapsePos, boolean collapseDinuc ) {
    public Set<String> getReadGroups() {
    public boolean isCollapseDinuc() {
    public boolean isCollapsePos() {
    public int getNReadGroups() {
    private RecalData getRecalData(String readGroup, int pos, int qual, char prevBase, char base) {
    public List<RecalData> getRecalData(String readGroup) {
    public int updateDataFromRead( String rg, SAMRecord read, int offset, char ref ) {
    */

    /**
     * The fasta, for comparison.
     */
    protected static IndexedFastaSequenceFile sequenceFile = null;

    CovariateCounter c;

    /**
     * Initialize the fasta.
     */
    @BeforeClass
    public static void initialize() throws FileNotFoundException {
        sequenceFile = new IndexedFastaSequenceFile( new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta") );
        GenomeLocParser.setupRefContigOrdering(sequenceFile);

    }

    @Before
    public void initializeBefore() {
        header = ArtificialSAMUtils.createArtificialSamHeader(2,0,247249719);
        readGroups.addAll(Arrays.asList(readGroup1, readGroup2));
        ArtificialSAMUtils.createDefaultReadGroup( header, readGroup1, "sample1" );
        ArtificialSAMUtils.createDefaultReadGroup( header, readGroup2, "sample2" ); 
        c = new CovariateCounter( readGroups, false, false, false );

        read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1, bases1, quals1);
        read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",1,1, bases2, quals2);
        read3 = ArtificialSAMUtils.createArtificialRead(header,"read3",1,1, bases1, quals3);
    }

    @Test
    public void testCovariateCounterSetup() {
        Assert.assertEquals("Number of read groups is wrong", c.getNReadGroups(), 2);
        Assert.assertEquals("Read group identities are wrong", c.getReadGroups(), readGroups);
        Assert.assertEquals("Incorrectly collapsed counter", c.isCollapseDinuc(), false);
        Assert.assertEquals("Incorrectly collapsed counter", c.isCollapsePos(), false);
    }

    @Test
    public void testOneRead() {
        for ( int i = 1; i < read1.getReadBases().length; i++ )
            c.updateDataFromRead(readGroup1, read1, i, (char)read1.getReadBases()[i]);
        c.printState();

        Assert.assertEquals("Incorrect mapping to recal bin", c.getRecalData(readGroup1, 0, quals1[0], 'A', (char)bases1[0]).N, 0);
        for ( int i = 1; i < bases1.length; i++ ) {
            RecalData datum = c.getRecalData(readGroup1, i, quals1[i], (char)bases1[i-1], (char)bases1[i]);
            System.out.printf("%s%n", datum);
            Assert.assertNotNull("Incorrect mapping to recal bin", datum);
            Assert.assertEquals("Bad mismatch count", datum.B, 0);
            Assert.assertEquals("Bad base count", datum.N, 1);
            Assert.assertEquals("Prevbase is bad", datum.dinuc.charAt(0), bases1[i-1]);
            Assert.assertEquals("Base is bad", datum.dinuc.charAt(1), bases1[i]);
            Assert.assertEquals("Qual is bad", datum.qual, quals1[i]);
        }
    }

    @Test
    public void testTwoReads() {
        for ( int i = 1; i < read1.getReadBases().length; i++ )
            c.updateDataFromRead(readGroup1, read1, i, (char)read1.getReadBases()[i]);
        for ( int i = 1; i < read2.getReadBases().length; i++ )
            c.updateDataFromRead(readGroup2, read2, i, (char)read2.getReadBases()[i]);
        c.printState();

        Assert.assertEquals("Incorrect mapping to recal bin", c.getRecalData(readGroup1, 0, quals1[0], 'A', (char)bases1[0]).N, 0);
        for ( int i = 1; i < bases1.length; i++ ) {
            RecalData datum = c.getRecalData(readGroup1, i, quals1[i], (char)bases1[i-1], (char)bases1[i]);
            System.out.printf("%s%n", datum);
            Assert.assertNotNull("Incorrect mapping to recal bin", datum);
            Assert.assertEquals("Bad mismatch count", datum.B, 0);
            Assert.assertEquals("Bad base count", datum.N, 1);
            Assert.assertEquals("Prevbase is bad", datum.dinuc.charAt(0), bases1[i-1]);
            Assert.assertEquals("Base is bad", datum.dinuc.charAt(1), bases1[i]);
            Assert.assertEquals("Qual is bad", datum.qual, quals1[i]);
        }
    }

    @Test
    public void testTwoReadsSameGroup() {
        for ( int i = 1; i < read1.getReadBases().length; i++ )
            c.updateDataFromRead(readGroup1, read1, i, (char)read1.getReadBases()[i]);
        for ( int i = 1; i < read2.getReadBases().length; i++ )
            c.updateDataFromRead(readGroup1, read1, i, (char)read1.getReadBases()[i]);
        c.printState();

        for ( int i = 1; i < bases1.length; i++ ) {
            RecalData datum = c.getRecalData(readGroup1, i, quals1[i], (char)bases1[i-1], (char)bases1[i]);
            System.out.printf("%s%n", datum);
            Assert.assertNotNull("Incorrect mapping to recal bin", datum);
            Assert.assertEquals("Bad mismatch count", datum.B, 0);
            Assert.assertEquals("Bad base count", datum.N, 2);
            Assert.assertEquals("Prevbase is bad", datum.dinuc.charAt(0), bases1[i-1]);
            Assert.assertEquals("Base is bad", datum.dinuc.charAt(1), bases1[i]);
            Assert.assertEquals("Qual is bad", datum.qual, quals1[i]);
        }
    }

    @Test
    public void testTwoReadsSameGroupNotIdentical() {
        for ( int i = 1; i < read1.getReadBases().length; i++ )
            c.updateDataFromRead(readGroup1, read1, i, (char)read1.getReadBases()[i]);
        for ( int i = 1; i < read3.getReadBases().length; i++ )
            c.updateDataFromRead(readGroup1, read3, i, (char)read3.getReadBases()[i]);
        c.printState();

        for ( int i = 1; i < bases1.length; i++ ) {
            RecalData datum = c.getRecalData(readGroup1, i, quals1[i], (char)bases1[i-1], (char)bases1[i]);
            System.out.printf("%s%n", datum);
            Assert.assertNotNull("Incorrect mapping to recal bin", datum);
            Assert.assertEquals("Bad mismatch count", datum.B, 0);
            Assert.assertEquals("Bad base count", datum.N, quals1[i] == quals3[i] ? 2 : 1);
            Assert.assertEquals("Prevbase is bad", datum.dinuc.charAt(0), bases1[i-1]);
            Assert.assertEquals("Base is bad", datum.dinuc.charAt(1), bases1[i]);
            Assert.assertEquals("Qual is bad", datum.qual, quals1[i]);
        }
    }

    @Test (expected = RuntimeException.class)
    public void testBadReadOffset() {
        byte bases[] = {'a', 't', 'c', 'g', 'a'};
        byte quals[] = {1, 2, 3, 4, 5};

        SAMRecord read = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1, bases, quals);

        c.updateDataFromRead(readGroup1, read, 0, (char)bases[0]);
    }
}