/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.haplotypeBAMWriter;

import net.sf.samtools.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class HaplotypeBAMWriterUnitTest extends BaseTest {
    private final static boolean DEBUG = false;
    final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);

    private GATKSAMRecord makeRead(final String baseString) {
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, 10);
        final byte[] bases = baseString.getBytes();
        read.setReadBases(bases.clone());
        read.setBaseQualities(Utils.dupBytes((byte)30, read.getReadLength()));
        return read;
    }

    private Haplotype makeHaplotype(final String bases, final String cigar) {
        final Haplotype hap = new Haplotype(bases.getBytes());
        hap.setCigar(TextCigarCodec.getSingleton().decode(cigar));
        return hap;
    }

    private static class MockBAMWriter implements SAMFileWriter {
        @Override
        public void addAlignment(SAMRecord alignment) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public SAMFileHeader getFileHeader() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public void close() {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    @Test
    public void testCreate() throws Exception {
        final SAMFileWriter writer = new MockBAMWriter();
        Assert.assertTrue(HaplotypeBAMWriter.create(HaplotypeBAMWriter.Type.CALLED_HAPLOTYPES, writer) instanceof CalledHaplotypeBAMWriter);
        Assert.assertTrue(HaplotypeBAMWriter.create(HaplotypeBAMWriter.Type.ALL_POSSIBLE_HAPLOTYPES, writer) instanceof AllHaplotypeBAMWriter);
    }


    //////////////////////////////////////////
    // Test HaplotypeBAMWriter.createReadAlignedToRef() //
    //////////////////////////////////////////

    @DataProvider(name = "ReadAlignedToRefData")
    public Object[][] makeReadAlignedToRefData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String hapBases = "ACTGAAGGTTCC";
        final Haplotype allM = makeHaplotype(hapBases, hapBases.length() + "M");

        // make sure we get back a cigar of the right length
        for ( int i = -1; i < hapBases.length(); i++ ) {
            final GATKSAMRecord read = makeRead(hapBases);
            if ( i != -1 ) read.getReadBases()[i] = (byte)'A';
            tests.add(new Object[]{read, allM, 10, 10, allM.getCigar().toString()});
        }

        // make sure insertions at the front are correctly handled
        for ( int padFront = 1; padFront < 10; padFront++ ) {
            final GATKSAMRecord read = makeRead(Utils.dupString("N", padFront) + hapBases);
            tests.add(new Object[]{read, allM, 10, 10, padFront + "I" + allM.getCigar().toString()});
        }

        // make sure insertions at the back are correctly handled
        for ( int padBack = 1; padBack < 10; padBack++ ) {
            final GATKSAMRecord read = makeRead(hapBases + Utils.dupString("N", padBack));
            tests.add(new Object[]{read, allM, 10, 10, allM.getCigar().toString() + padBack + "I"});
        }

        // make sure refStart and hapStart are respected
        for ( int refStart = 1; refStart < 10; refStart++ ) {
            for ( int hapStart = refStart; hapStart < 10 + refStart; hapStart++ ) {
                final Haplotype hap = new Haplotype(allM.getBases());
                hap.setCigar(allM.getCigar());
                hap.setAlignmentStartHapwrtRef(hapStart);

                final GATKSAMRecord read = makeRead(new String(hap.getBases()));
                tests.add(new Object[]{read, hap, refStart, refStart + hapStart, allM.getCigar().toString()});
            }
        }

        // example case of bad alignment because SW doesn't necessarily left-align indels
        {
            final String hap = "ACTGTGGGTTCCTCTTATTTTATTTCTACATCAATGTTCATATTTAACTTATTATTTTATCTTATTTTTAAATTTCTTTTATGTTGAGCCTTGATGAAAGCCATAGGTTCTCTCATATAATTGTATGTGTATGTATGTATATGTACATAATATATACATATATGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATATATACATATATG";
            final String hapCigar = "399M";
            final String readBases = "ATGTACATAATATATACATATATGTATATGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATAT";
            final GATKSAMRecord read = makeRead(readBases);
            final int refStart = 10130100;
            final int hapStart = 500;
            final String badCigar = "31M6D211M";
            final String goodCigar = "28M6D214M";
            final Haplotype badHap = new Haplotype(hap.getBytes());
            badHap.setCigar(TextCigarCodec.getSingleton().decode(hapCigar));
            badHap.setAlignmentStartHapwrtRef(hapStart);

            final int expectedPos = 10130740;
            tests.add(new Object[]{read, badHap, refStart, expectedPos, goodCigar});
        }

        return tests.toArray(new Object[][]{});
    }



    @Test(dataProvider = "ReadAlignedToRefData", enabled = true)
    public void testReadAlignedToRef(final GATKSAMRecord read, final Haplotype haplotype, final int refStart, final int expectedReadStart, final String expectedReadCigar) throws Exception {
        final HaplotypeBAMWriter writer = new CalledHaplotypeBAMWriter(new MockBAMWriter());
        final GATKSAMRecord originalReadCopy = (GATKSAMRecord)read.clone();

        if ( expectedReadCigar == null ) {
            Assert.assertNull(writer.createReadAlignedToRef(read, haplotype, refStart));
        } else {
            final Cigar expectedCigar = TextCigarCodec.getSingleton().decode(expectedReadCigar);
            final GATKSAMRecord alignedRead = writer.createReadAlignedToRef(read, haplotype, refStart);

            Assert.assertEquals(alignedRead.getReadName(), originalReadCopy.getReadName());
            Assert.assertEquals(alignedRead.getAlignmentStart(), expectedReadStart);
            Assert.assertEquals(alignedRead.getReadBases(), originalReadCopy.getReadBases());
            Assert.assertEquals(alignedRead.getBaseQualities(), originalReadCopy.getBaseQualities());
            Assert.assertEquals(alignedRead.getAlignmentStart(), expectedReadStart);
            Assert.assertEquals(alignedRead.getCigar(), expectedCigar);
            Assert.assertNotNull(alignedRead.getAttribute("HC"));
        }

        Assert.assertEquals(read, originalReadCopy, "createReadAlignedToRef seems be modifying the original read!");
    }

    private static class Mutation implements Comparable<Mutation> {
        int pos, len;
        CigarOperator operator;

        private Mutation(int pos, int len, CigarOperator operator) {
            this.pos = pos;
            this.len = len;
            this.operator = operator;
        }
        public int getNMismatches() { return len; }

        @Override
        public int compareTo(Mutation o) {
            return Integer.valueOf(pos).compareTo(o.pos);
        }

        private String apply(final String seq) {
            switch ( operator ) {
                case M:
                    final byte[] bases = seq.getBytes();
                    if ( pos < seq.length() )
                        bases[pos] = (byte)(bases[pos] == 'A' ? 'C' : 'A');
                    return new String(bases);
                case I: {
                    final String prefix = seq.substring(0, pos);
                    final String postfix = seq.substring(pos, seq.length());
                    return prefix + "GTCAGTTA".substring(0, len) + postfix;
                } case D: {
                    final String prefix = seq.substring(0, pos);
                    final String postfix = seq.substring(pos + len, seq.length());
                    return prefix + postfix;
                }default:
                    throw new IllegalStateException("Unexpected operator " + operator);
            }
        }
    }

    private static class MutatedSequence {
        int numMismatches;
        String seq;

        private MutatedSequence(int numMismatches, String seq) {
            this.numMismatches = numMismatches;
            this.seq = seq;
        }
    }

    private MutatedSequence mutateSequence(final String hapIn, final List<Mutation> mutations) {
        Collections.sort(mutations);
        int mismatches = 0;
        String hap = hapIn;
        for ( final Mutation mut : mutations ) {
            hap = mut.apply(hap);
            mismatches += mut.getNMismatches();
        }
        return new MutatedSequence(mismatches, hap);
    }

    @DataProvider(name = "ComplexReadAlignedToRef")
    public Object[][] makeComplexReadAlignedToRef() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<Mutation> allMutations = Arrays.asList(
                new Mutation(1, 1, CigarOperator.M),
                new Mutation(2, 1, CigarOperator.M),
                new Mutation(3, 1, CigarOperator.I),
                new Mutation(7, 1, CigarOperator.D)
        );

        int i = 0;
        final String referenceBases  = "ACTGACTGACTG";
        final String paddedReference = "NNNN" + referenceBases + "NNNN";
        for ( final List<Mutation> mutations : Utils.makePermutations(allMutations, 3, false) ) {
            final MutatedSequence hap = mutateSequence(referenceBases, mutations);
            final Haplotype haplotype = new Haplotype(hap.seq.getBytes());
            final SWPairwiseAlignment align = new SWPairwiseAlignment(paddedReference.getBytes(), hap.seq.getBytes());
            haplotype.setAlignmentStartHapwrtRef(align.getAlignmentStart2wrt1());
            haplotype.setCigar(align.getCigar());

            for ( final List<Mutation> readMutations : Utils.makePermutations(allMutations, 3, false) ) {
                final MutatedSequence readBases = mutateSequence(hap.seq, readMutations);
                final GATKSAMRecord read = makeRead(readBases.seq);
                tests.add(new Object[]{i++, read, paddedReference, haplotype, hap.numMismatches + readBases.numMismatches});
            }
        }

        // for convenient testing of a single failing case
        //tests.add(new Object[]{makeRead("ACCGGGACTGACTG"), reference, makeHaplotype("AAAGGACTGACTG", "1M1I11M"), 2});

        return tests.toArray(new Object[][]{});
    }


    @Test(dataProvider = "ComplexReadAlignedToRef", enabled = !DEBUG)
    public void testReadAlignedToRefComplexAlignment(final int testIndex, final GATKSAMRecord read, final String reference, final Haplotype haplotype, final int expectedMaxMismatches) throws Exception {
        final HaplotypeBAMWriter writer = new CalledHaplotypeBAMWriter(new MockBAMWriter());
        final GATKSAMRecord alignedRead = writer.createReadAlignedToRef(read, haplotype, 1);
        if ( alignedRead != null ) {
            final int mismatches = AlignmentUtils.getMismatchCount(alignedRead, reference.getBytes(), alignedRead.getAlignmentStart() - 1).numMismatches;
            Assert.assertTrue(mismatches <= expectedMaxMismatches,
                    "Alignment of read to ref looks broken.  Expected at most " + expectedMaxMismatches + " but saw " + mismatches
                            + " for readBases " + new String(read.getReadBases()) + " with cigar " + read.getCigar() + " reference " + reference + " haplotype "
                            + haplotype + " with cigar " + haplotype.getCigar() + " aligned read cigar " + alignedRead.getCigarString() + " @ " + alignedRead.getAlignmentStart());
        }
    }
}
