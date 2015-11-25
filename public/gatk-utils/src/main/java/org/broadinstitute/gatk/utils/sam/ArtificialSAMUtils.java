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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;

import java.io.File;
import java.util.*;

/**
 * @author aaron
 * @version 1.0
 */
public class ArtificialSAMUtils {
    public static final int DEFAULT_READ_LENGTH = 50;

    /**
     * create an artificial sam file
     *
     * @param filename            the filename to write to
     * @param numberOfChromosomes the number of chromosomes
     * @param startingChromosome  where to start counting
     * @param chromosomeSize      how large each chromosome is
     * @param readsPerChomosome   how many reads to make in each chromosome.  They'll be aligned from position 1 to x (which is the number of reads)
     */
    public static void createArtificialBamFile(String filename, int numberOfChromosomes, int startingChromosome, int chromosomeSize, int readsPerChomosome) {
        SAMFileHeader header = createArtificialSamHeader(numberOfChromosomes, startingChromosome, chromosomeSize);
        File outFile = new File(filename);

        SAMFileWriter out = new SAMFileWriterFactory().makeBAMWriter(header, true, outFile);

        for (int x = startingChromosome; x < startingChromosome + numberOfChromosomes; x++) {
            for (int readNumber = 1; readNumber < readsPerChomosome; readNumber++) {
                out.addAlignment(createArtificialRead(header, "Read_" + readNumber, x - startingChromosome, readNumber, DEFAULT_READ_LENGTH));
            }
        }

        out.close();
    }

    /**
     * create an artificial sam file
     *
     * @param filename            the filename to write to
     * @param numberOfChromosomes the number of chromosomes
     * @param startingChromosome  where to start counting
     * @param chromosomeSize      how large each chromosome is
     * @param readsPerChomosome   how many reads to make in each chromosome.  They'll be aligned from position 1 to x (which is the number of reads)
     */
    public static void createArtificialSamFile(String filename, int numberOfChromosomes, int startingChromosome, int chromosomeSize, int readsPerChomosome) {
        SAMFileHeader header = createArtificialSamHeader(numberOfChromosomes, startingChromosome, chromosomeSize);
        File outFile = new File(filename);

        SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header, false, outFile);

        for (int x = startingChromosome; x < startingChromosome + numberOfChromosomes; x++) {
            for (int readNumber = 1; readNumber <= readsPerChomosome; readNumber++) {
                out.addAlignment(createArtificialRead(header, "Read_" + readNumber, x - startingChromosome, readNumber, 100));
            }
        }

        out.close();
    }

    /**
     * Creates an artificial sam header, matching the parameters, chromosomes which will be labeled chr1, chr2, etc
     *
     * @param numberOfChromosomes the number of chromosomes to create
     * @param startingChromosome  the starting number for the chromosome (most likely set to 1)
     * @param chromosomeSize      the length of each chromosome
     * @return
     */
    public static SAMFileHeader createArtificialSamHeader(int numberOfChromosomes, int startingChromosome, int chromosomeSize) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        // make up some sequence records
        for (int x = startingChromosome; x < startingChromosome + numberOfChromosomes; x++) {
            SAMSequenceRecord rec = new SAMSequenceRecord("chr" + (x), chromosomeSize /* size */);
            rec.setSequenceLength(chromosomeSize);
            dict.addSequence(rec);
        }
        header.setSequenceDictionary(dict);
        return header;
    }

    /**
     * Creates an artificial sam header based on the sequence dictionary dict
     *
     * @return a new sam header
     */
    public static SAMFileHeader createArtificialSamHeader(final SAMSequenceDictionary dict) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);
        header.setSequenceDictionary(dict);
        return header;
    }

    /**
     * Creates an artificial sam header with standard test parameters
     *
     * @return the sam header
     */
    public static SAMFileHeader createArtificialSamHeader() {
        return createArtificialSamHeader(1, 1, 1000000);
    }

    /**
     * setup a default read group for a SAMFileHeader
     *
     * @param header      the header to set
     * @param readGroupID the read group ID tag
     * @param sampleName  the sample name
     * @return the adjusted SAMFileHeader
     */
    public static SAMFileHeader createDefaultReadGroup(SAMFileHeader header, String readGroupID, String sampleName) {
        SAMReadGroupRecord rec = new SAMReadGroupRecord(readGroupID);
        rec.setSample(sampleName);
        List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        readGroups.add(rec);
        header.setReadGroups(readGroups);
        return header;
    }

    /**
     * setup read groups for the specified read groups and sample names
     *
     * @param header       the header to set
     * @param readGroupIDs the read group ID tags
     * @param sampleNames  the sample names
     * @return the adjusted SAMFileHeader
     */
    public static SAMFileHeader createEnumeratedReadGroups(SAMFileHeader header, List<String> readGroupIDs, List<String> sampleNames) {
        if (readGroupIDs.size() != sampleNames.size()) {
            throw new ReviewedGATKException("read group count and sample name count must be the same");
        }

        List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();

        int x = 0;
        for (; x < readGroupIDs.size(); x++) {
            SAMReadGroupRecord rec = new SAMReadGroupRecord(readGroupIDs.get(x));
            rec.setSample(sampleNames.get(x));
            readGroups.add(rec);
        }
        header.setReadGroups(readGroups);
        return header;
    }


    /**
     * Create an artificial read based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param length         the length of the read
     * @return the artificial read
     */
    public static GATKSAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        if ((refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart != SAMRecord.NO_ALIGNMENT_START) ||
                (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart == SAMRecord.NO_ALIGNMENT_START))
            throw new ReviewedGATKException("Invalid alignment start for artificial read, start = " + alignmentStart);
        GATKSAMRecord record = new GATKSAMRecord(header);
        record.setReadName(name);
        record.setReferenceIndex(refIndex);
        record.setAlignmentStart(alignmentStart);
        List<CigarElement> elements = new ArrayList<CigarElement>();
        elements.add(new CigarElement(length, CigarOperator.characterToEnum('M')));
        record.setCigar(new Cigar(elements));
        record.setProperPairFlag(false);

        // our reads and quals are all 'A's by default
        byte[] c = new byte[length];
        byte[] q = new byte[length];
        for (int x = 0; x < length; x++)
            c[x] = q[x] = 'A';
        record.setReadBases(c);
        record.setBaseQualities(q);

        if (refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            record.setReadUnmappedFlag(true);
        }

        return record;
    }

    /**
     * Create an artificial read based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @return the artificial read
     */
    public static GATKSAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual) {
        if (bases.length != qual.length) {
            throw new ReviewedGATKException("Passed in read string is different length then the quality array");
        }
        GATKSAMRecord rec = createArtificialRead(header, name, refIndex, alignmentStart, bases.length);
        rec.setReadBases(bases);
        rec.setBaseQualities(qual);
        rec.setReadGroup(new GATKSAMReadGroupRecord("x"));
        if (refIndex == -1) {
            rec.setReadUnmappedFlag(true);
        }

        return rec;
    }

    /**
     * Create an artificial read based on the parameters
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @param cigar          the cigar string of the read
     * @return the artificial read
     */
    public static GATKSAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual, String cigar) {
        GATKSAMRecord rec = createArtificialRead(header, name, refIndex, alignmentStart, bases, qual);
        rec.setCigarString(cigar);
        return rec;
    }

    /**
     * Create an artificial read with the following default parameters :
     * header:
     * numberOfChromosomes = 1
     * startingChromosome = 1
     * chromosomeSize = 1000000
     * read:
     * name = "default_read"
     * refIndex = 0
     * alignmentStart = 1
     *
     * @param bases the sequence of the read
     * @param qual  the qualities of the read
     * @param cigar the cigar string of the read
     * @return the artificial read
     */
    public static GATKSAMRecord createArtificialRead(byte[] bases, byte[] qual, String cigar) {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();
        return ArtificialSAMUtils.createArtificialRead(header, "default_read", 0, 10000, bases, qual, cigar);
    }

    public static GATKSAMRecord createArtificialRead(Cigar cigar) {
        int length = cigar.getReadLength();
        byte [] base = {'A'};
        byte [] qual = {30};
        byte [] bases = Utils.arrayFromArrayWithLength(base, length);
        byte [] quals = Utils.arrayFromArrayWithLength(qual, length);
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();
        return ArtificialSAMUtils.createArtificialRead(header, "default_read", 0, 10000, bases, quals, cigar.toString());
    }

    
    public final static List<GATKSAMRecord> createPair(SAMFileHeader header, String name, int readLen, int leftStart, int rightStart, boolean leftIsFirst, boolean leftIsNegative) {
        GATKSAMRecord left = ArtificialSAMUtils.createArtificialRead(header, name, 0, leftStart, readLen);
        GATKSAMRecord right = ArtificialSAMUtils.createArtificialRead(header, name, 0, rightStart, readLen);

        left.setReadPairedFlag(true);
        right.setReadPairedFlag(true);

        left.setProperPairFlag(true);
        right.setProperPairFlag(true);

        left.setFirstOfPairFlag(leftIsFirst);
        right.setFirstOfPairFlag(!leftIsFirst);

        left.setReadNegativeStrandFlag(leftIsNegative);
        left.setMateNegativeStrandFlag(!leftIsNegative);
        right.setReadNegativeStrandFlag(!leftIsNegative);
        right.setMateNegativeStrandFlag(leftIsNegative);

        left.setMateAlignmentStart(right.getAlignmentStart());
        right.setMateAlignmentStart(left.getAlignmentStart());

        left.setMateReferenceIndex(0);
        right.setMateReferenceIndex(0);

        int isize = rightStart + readLen - leftStart;
        left.setInferredInsertSize(isize);
        right.setInferredInsertSize(-isize);

        return Arrays.asList(left, right);
    }

    /**
     * Create a collection of identical artificial reads based on the parameters.  The cigar string for each
     * read will be *M, where * is the length of the read.
     *
     * Useful for testing things like positional downsampling where you care only about the position and
     * number of reads, and not the other attributes.
     *
     * @param stackSize      number of identical reads to create
     * @param header         the SAM header to associate each read with
     * @param name           name associated with each read
     * @param refIndex       the reference index, i.e. what chromosome to associate them with
     * @param alignmentStart where to start each alignment
     * @param length         the length of each read
     *
     * @return a collection of stackSize reads all sharing the above properties
     */
    public static Collection<GATKSAMRecord> createStackOfIdenticalArtificialReads( int stackSize, SAMFileHeader header, String name, int refIndex, int alignmentStart, int length ) {
        Collection<GATKSAMRecord> stack = new ArrayList<GATKSAMRecord>(stackSize);
        for ( int i = 1; i <= stackSize; i++ ) {
            stack.add(createArtificialRead(header, name, refIndex, alignmentStart, length));
        }
        return stack;
    }

    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr the chromosome (reference ID) to start from
     * @param endingChr   the id to end with
     * @param readCount   the number of reads per chromosome
     * @return GATKSAMIterator representing the specified amount of fake data
     */
    public static GATKSAMIterator mappedReadIterator(int startingChr, int endingChr, int readCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, 0, header);
    }

    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr       the chromosome (reference ID) to start from
     * @param endingChr         the id to end with
     * @param readCount         the number of reads per chromosome
     * @param unmappedReadCount the count of unmapped reads to place at the end of the iterator, like in a sorted bam file
     * @return GATKSAMIterator representing the specified amount of fake data
     */
    public static GATKSAMIterator mappedAndUnmappedReadIterator(int startingChr, int endingChr, int readCount, int unmappedReadCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, unmappedReadCount, header);
    }

    /**
     * create an ArtificialSAMQueryIterator containing the specified read piles
     *
     * @param startingChr the chromosome (reference ID) to start from
     * @param endingChr   the id to end with
     * @param readCount   the number of reads per chromosome
     * @return GATKSAMIterator representing the specified amount of fake data
     */
    public static ArtificialSAMQueryIterator queryReadIterator(int startingChr, int endingChr, int readCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, 0, header);
    }

    /**
     * create an ArtificialSAMQueryIterator containing the specified read piles
     *
     * @param startingChr       the chromosome (reference ID) to start from
     * @param endingChr         the id to end with
     * @param readCount         the number of reads per chromosome
     * @param unmappedReadCount the count of unmapped reads to place at the end of the iterator, like in a sorted bam file
     * @return GATKSAMIterator representing the specified amount of fake data
     */
    public static GATKSAMIterator queryReadIterator(int startingChr, int endingChr, int readCount, int unmappedReadCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, unmappedReadCount, header);
    }

    /**
     * Create an iterator containing the specified reads
     *
     * @param reads the reads
     * @return iterator for the reads
     */
    public static GATKSAMIterator createReadIterator(SAMRecord... reads) {
        return createReadIterator(Arrays.asList(reads));
    }

    /**
     * Create an iterator containing the specified reads
     *
     * @param reads the reads
     * @return iterator for the reads
     */
    public static GATKSAMIterator createReadIterator(List<SAMRecord> reads) {
        final Iterator<SAMRecord> iter = reads.iterator();
        return new GATKSAMIterator() {
            @Override public void close() {}
            @Override public Iterator<SAMRecord> iterator() { return iter; }
            @Override public boolean hasNext() { return iter.hasNext(); }
            @Override public SAMRecord next() { return iter.next(); }
            @Override public void remove() { iter.remove(); }
        };
    }

    private final static int ranIntInclusive(Random ran, int start, int stop) {
        final int range = stop - start;
        return ran.nextInt(range) + start;
    }

    /**
     * Creates a read backed pileup containing up to pileupSize reads at refID 0 from header at loc with
     * reads created that have readLen bases.  Pairs are sampled from a gaussian distribution with mean insert
     * size of insertSize and variation of insertSize / 10.  The first read will be in the pileup, and the second
     * may be, depending on where this sampled insertSize puts it.
     *
     * @param header
     * @param loc
     * @param readLen
     * @param insertSize
     * @param pileupSize
     * @return
     */
    public static ReadBackedPileup createReadBackedPileup(final SAMFileHeader header, final GenomeLoc loc, final int readLen, final int insertSize, final int pileupSize) {
        final Random ran = new Random();
        final boolean leftIsFirst = true;
        final boolean leftIsNegative = false;
        final int insertSizeVariation = insertSize / 10;
        final int pos = loc.getStart();

        final List<PileupElement> pileupElements = new ArrayList<PileupElement>();
        for (int i = 0; i < pileupSize / 2; i++) {
            final String readName = "read" + i;
            final int leftStart = ranIntInclusive(ran, 1, pos);
            final int fragmentSize = (int) (ran.nextGaussian() * insertSizeVariation + insertSize);
            final int rightStart = leftStart + fragmentSize - readLen;

            if (rightStart <= 0) continue;

            List<GATKSAMRecord> pair = createPair(header, readName, readLen, leftStart, rightStart, leftIsFirst, leftIsNegative);
            final GATKSAMRecord left = pair.get(0);
            final GATKSAMRecord right = pair.get(1);

            pileupElements.add(LocusIteratorByState.createPileupForReadAndOffset(left, pos - leftStart));

            if (pos >= right.getAlignmentStart() && pos <= right.getAlignmentEnd()) {
                pileupElements.add(LocusIteratorByState.createPileupForReadAndOffset(right, pos - rightStart));
            }
        }

        Collections.sort(pileupElements);
        return new ReadBackedPileupImpl(loc, pileupElements);
    }
}
