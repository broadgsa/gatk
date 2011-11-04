package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;

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
    public static void createArtificialBamFile( String filename, int numberOfChromosomes, int startingChromosome, int chromosomeSize, int readsPerChomosome ) {
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
    public static void createArtificialSamFile( String filename, int numberOfChromosomes, int startingChromosome, int chromosomeSize, int readsPerChomosome ) {
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
     *
     * @return
     */
    public static SAMFileHeader createArtificialSamHeader( int numberOfChromosomes, int startingChromosome, int chromosomeSize ) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(net.sf.samtools.SAMFileHeader.SortOrder.coordinate);
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        // make up some sequence records
        for (int x = startingChromosome; x < startingChromosome + numberOfChromosomes; x++) {
            SAMSequenceRecord rec = new SAMSequenceRecord("chr" + ( x ), chromosomeSize /* size */);
            rec.setSequenceLength(chromosomeSize);
            dict.addSequence(rec);
        }
        header.setSequenceDictionary(dict);
        return header;
    }

    /**
     * setup a default read group for a SAMFileHeader
     *
     * @param header      the header to set
     * @param readGroupID the read group ID tag
     * @param sampleName  the sample name
     *
     * @return the adjusted SAMFileHeader
     */
    public static SAMFileHeader createDefaultReadGroup( SAMFileHeader header, String readGroupID, String sampleName ) {
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
     *
     * @return the adjusted SAMFileHeader
     */
    public static SAMFileHeader createEnumeratedReadGroups( SAMFileHeader header, List<String> readGroupIDs, List<String> sampleNames ) {
        if (readGroupIDs.size() != sampleNames.size()) {
            throw new ReviewedStingException("read group count and sample name count must be the same");
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
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param length         the length of the read
     *
     * @return the artificial read
     */
    public static GATKSAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        if( (refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart != SAMRecord.NO_ALIGNMENT_START) ||
                (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart == SAMRecord.NO_ALIGNMENT_START) )
            throw new ReviewedStingException("Invalid alignment start for artificial read, start = " + alignmentStart);
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
     *
     * @return the artificial read
     */
    public static GATKSAMRecord createArtificialRead( SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual ) {
        if (bases.length != qual.length) {
            throw new ReviewedStingException("Passed in read string is different length then the quality array");
        }
        GATKSAMRecord rec = createArtificialRead(header, name, refIndex, alignmentStart, bases.length);
        rec.setReadBases(bases);
        rec.setBaseQualities(qual);
        if (refIndex == -1) {
            rec.setReadUnmappedFlag(true);
        }

        return rec;
    }

    public final static List<GATKSAMRecord> createPair(SAMFileHeader header, String name, int readLen, int leftStart, int rightStart, boolean leftIsFirst, boolean leftIsNegative) {
        GATKSAMRecord left = ArtificialSAMUtils.createArtificialRead(header, name, 0, leftStart, readLen);
        GATKSAMRecord right = ArtificialSAMUtils.createArtificialRead(header, name, 0, rightStart, readLen);

        left.setReadPairedFlag(true);
        right.setReadPairedFlag(true);

        left.setProperPairFlag(true);
        right.setProperPairFlag(true);

        left.setFirstOfPairFlag(leftIsFirst);
        right.setFirstOfPairFlag(! leftIsFirst);

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
     * create an iterator containing the specified read piles
     *
     * @param startingChr the chromosome (reference ID) to start from
     * @param endingChr   the id to end with
     * @param readCount   the number of reads per chromosome
     *
     * @return StingSAMIterator representing the specified amount of fake data
     */
    public static StingSAMIterator mappedReadIterator( int startingChr, int endingChr, int readCount ) {
        SAMFileHeader header = createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, 0, header);
    }

    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr       the chromosome (reference ID) to start from
     * @param endingChr         the id to end with
     * @param readCount         the number of reads per chromosome
     * @param unmappedReadCount the count of unmapped reads to place at the end of the iterator, like in a sorted bam file
     *
     * @return StingSAMIterator representing the specified amount of fake data
     */
    public static StingSAMIterator mappedAndUnmappedReadIterator( int startingChr, int endingChr, int readCount, int unmappedReadCount ) {
        SAMFileHeader header = createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, unmappedReadCount, header);
    }

    /**
     * create an ArtificialSAMQueryIterator containing the specified read piles
     *
     * @param startingChr the chromosome (reference ID) to start from
     * @param endingChr   the id to end with
     * @param readCount   the number of reads per chromosome
     *
     * @return StingSAMIterator representing the specified amount of fake data
     */
    public static ArtificialSAMQueryIterator queryReadIterator( int startingChr, int endingChr, int readCount ) {
        SAMFileHeader header = createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, 0, header);
    }

    /**
     * create an ArtificialSAMQueryIterator containing the specified read piles
     *
     * @param startingChr       the chromosome (reference ID) to start from
     * @param endingChr         the id to end with
     * @param readCount         the number of reads per chromosome
     * @param unmappedReadCount the count of unmapped reads to place at the end of the iterator, like in a sorted bam file
     *
     * @return StingSAMIterator representing the specified amount of fake data
     */
    public static StingSAMIterator queryReadIterator( int startingChr, int endingChr, int readCount, int unmappedReadCount ) {
        SAMFileHeader header = createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, unmappedReadCount, header);
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
        for ( int i = 0; i < pileupSize / 2; i++ ) {
            final String readName = "read" + i;
            final int leftStart = ranIntInclusive(ran, 1, pos);
            final int fragmentSize = (int)(ran.nextGaussian() * insertSizeVariation + insertSize);
            final int rightStart = leftStart + fragmentSize - readLen;

            if ( rightStart <= 0 ) continue;

            List<GATKSAMRecord> pair = createPair(header, readName, readLen, leftStart, rightStart, leftIsFirst, leftIsNegative);
            final GATKSAMRecord left = pair.get(0);
            final GATKSAMRecord right = pair.get(1);

            pileupElements.add(new PileupElement(left, pos - leftStart));

            if ( pos >= right.getAlignmentStart() && pos <= right.getAlignmentEnd() ) {
                pileupElements.add(new PileupElement(right, pos - rightStart));
            }
        }

        Collections.sort(pileupElements);
        return new ReadBackedPileupImpl(loc, pileupElements);
    }
}
