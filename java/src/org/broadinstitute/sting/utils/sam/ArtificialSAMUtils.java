package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.iterators.QueryIterator;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

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
            throw new StingException("read group count and sample name count must be the same");
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
     *
     * @return the artificial read
     */
    public static SAMRecord createArtificialRead( SAMFileHeader header, String name, int refIndex, int alignmentStart, int length ) {
        if( (refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart != SAMRecord.NO_ALIGNMENT_START) ||
            (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart == SAMRecord.NO_ALIGNMENT_START) )
            throw new StingException("Invalid alignment start for artificial read, start = " + alignmentStart);
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceIndex(refIndex);
        record.setAlignmentStart(alignmentStart);
        List<CigarElement> elements = new ArrayList<CigarElement>();
        elements.add(new CigarElement(length, CigarOperator.characterToEnum('M')));
        record.setCigar(new Cigar(elements));
        record.setProperPairFlag(false);
        // TODO: add to code after checking why it breaks the clean read injector test
        /*byte[] c = new byte[length];
        for (int x = 0; x < length; x++)
            c[x] = 'A';  */
        //record.setReadBases(c);
        if (refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            record.setReadUmappedFlag(true);
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
    public static SAMRecord createArtificialRead( SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual ) {
        if (bases.length != qual.length) {
            throw new StingException("Passed in read string is different length then the quality array");
        }
        SAMRecord rec = createArtificialRead(header, name, refIndex, alignmentStart, bases.length);
        rec.setReadBases(bases);
        rec.setBaseQualities(qual);
        if (refIndex == -1) {
            rec.setReadUmappedFlag(true);
        }
        return rec;
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
    public static QueryIterator mappedReadIterator( int startingChr, int endingChr, int readCount ) {
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
    public static ArtificialSAMIterator mappedAndUnmappedReadIterator( int startingChr, int endingChr, int readCount, int unmappedReadCount ) {
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
    public static ArtificialSAMQueryIterator queryReadIterator( int startingChr, int endingChr, int readCount, int unmappedReadCount ) {
        SAMFileHeader header = createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, unmappedReadCount, header);
    }
}
