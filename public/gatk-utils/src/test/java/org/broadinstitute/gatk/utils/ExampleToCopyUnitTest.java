/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils;


// the imports for unit testing.


import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.ArtificialBAMBuilder;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class ExampleToCopyUnitTest extends BaseTest {
    // example genome loc parser for this test, can be deleted if you don't use the reference
    private GenomeLocParser genomeLocParser;

    // example fasta index file, can be deleted if you don't use the reference
    private IndexedFastaSequenceFile seq;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    /**
     * Combinatorial unit test data provider example.
     *
     * Creates data for testMyData test function, containing two arguments, start and size at each value
     *
     * @return Object[][] for testng DataProvider
     */
    @DataProvider(name = "MyDataProvider")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( final int start : Arrays.asList(1, 10, 100) ) {
            for ( final int size : Arrays.asList(1, 10, 100, 1000) ) {
                tests.add(new Object[]{start, size});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "MyDataProvider")
    public void testMyData(final int start, final int size) {
        // adaptor this code to do whatever testing you want given the arguments start and size
        Assert.assertTrue(start >= 0);
        Assert.assertTrue(size >= 0);
    }

    /**
     * DataProvider example using a class-based data structure
     */
    private class MyDataProviderClass extends TestDataProvider {
        private int start;
        private int size;

        private MyDataProviderClass(int start, int size) {
            super(MyDataProviderClass.class);
            this.start = start;
            this.size = size;
        }
    }

    @DataProvider(name = "MyClassBasedDataProvider")
    public Object[][] makeMyDataProviderClass() {
        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( final int start : Arrays.asList(1, 10, 100) ) {
            for ( final int size : Arrays.asList(1, 10, 100, 1000) ) {
                new MyDataProviderClass(start, size);
            }
        }

        return TestDataProvider.getTests(MyDataProviderClass.class);
    }

    /**
     * Example testng test using MyClassBasedDataProvider
     */
    @Test(dataProvider = "MyClassBasedDataProvider")
    public void testMyDataProviderClass(MyDataProviderClass testSpec) {
        // adaptor this code to do whatever testing you want given the arguments start and size
        Assert.assertTrue(testSpec.start >= 0);
        Assert.assertTrue(testSpec.size >= 0);
    }

    /**
     * A unit test that creates an artificial read for testing some code that uses reads
     */
    @Test()
    public void testWithARead() {
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, 10);
        Assert.assertEquals(read.getReadLength(), 10);
        // TODO -- add some tests here using read
    }

    /**
     * A unit test that creates a GenomeLoc for testing
     */
    @Test()
    public void testWithAGenomeLoc() {
        final GenomeLoc loc = genomeLocParser.createGenomeLoc("1", 1, 10);
        Assert.assertEquals(loc.size(), 10);
        // TODO -- add some tests here using the loc
    }

    /**
     * A unit test that creates an artificial read for testing some code that uses reads
     *
     * Note that effective creation of RBPs isn't so good.  If you need pileups of specific properties, you shoud
     * look into building them yourself as in the example below
     */
    @Test()
    public void testWithAPileup() {
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final GenomeLoc myLocation = genomeLocParser.createGenomeLoc("1", 10);
        final ReadBackedPileup pileup = ArtificialSAMUtils.createReadBackedPileup(header, myLocation, 10, 400, 10);
        Assert.assertFalse(pileup.isEmpty());
        // TODO -- add some tests here using pileup
    }

    /**
     * A unit test that creates an artificial read for testing some code that uses reads
     *
     * Builds the pileup from scratch to have specific properties
     */
    @Test()
    public void testBuildingAPileupWithSpecificProperties() {
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final GenomeLoc myLocation = genomeLocParser.createGenomeLoc("1", 10);

        final int pileupSize = 100;
        final int readLength = 10;
        final List<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        for ( int i = 0; i < pileupSize; i++ ) {
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte)'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other base is a C

            // set the read's bases and quals
            read.setReadBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            reads.add(read);
        }

        // create a pileup with all reads having offset 0
        final ReadBackedPileup pileup = new ReadBackedPileupImpl(myLocation, reads, 0);
        // TODO -- add some tests here using pileup

        // this code ensures that the pileup example is correct.  Can be deleted
        Assert.assertEquals(pileup.getNumberOfElements(), pileupSize);
        int nA = 0, nC = 0;
        for ( final PileupElement p : pileup ) {
            if ( p.getBase() == 'A' ) nA++;
            if ( p.getBase() == 'C' ) nC++;
        }
        Assert.assertEquals(nA, pileupSize / 2);
        Assert.assertEquals(nC, pileupSize / 2);

    }

    /**
     * A unit test that creates an artificial read for testing some code that uses reads
     */
    @Test()
    public void testWithBAMFile() {
        // create a fake BAM file, and iterate through it
        final ArtificialBAMBuilder bamBuilder = new ArtificialBAMBuilder(seq, 20, 10);
        final File bam = bamBuilder.makeTemporarilyBAMFile();
        final SAMFileReader reader = new SAMFileReader(bam);

        final Iterator<SAMRecord> bamIt = reader.iterator();
        while ( bamIt.hasNext() ) {
            final SAMRecord read = bamIt.next(); // all reads are actually GATKSAMRecords
            // TODO -- add some tests that use reads from a BAM
        }
    }

    /**
     * Test code that creates VariantContexts
     */
    @Test()
    public void testWithVariantContext() throws Exception {
        final List<Allele> alleles = Arrays.asList(Allele.create("A", true), Allele.create("C"));
        final VariantContext vc = new VariantContextBuilder("test", "1", 10, 10, alleles).make();
        Assert.assertTrue(vc.getAlleles().size() >= 0);
        // TODO -- add some tests that use VariantContext
    }
}