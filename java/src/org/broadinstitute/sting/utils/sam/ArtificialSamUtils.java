package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.*;

import java.io.File;
import java.util.*;

import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.Reads;

/**
 *
 * User: aaron
 * Date: May 21, 2009
 * Time: 2:57:48 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date May 21, 2009
 * <p/>
 * Class samUtils
 * <p/>
 * A collection of utilities for making sam and bam files.  Mostly these are for creating sam and bam files for testing purposes.
 */
public class ArtificialSamUtils {
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

        SAMFileWriter out = new SAMFileWriterFactory().makeBAMWriter(header, false, outFile);

        for (int x = startingChromosome; x < startingChromosome + numberOfChromosomes; x++) {
            for (int readNumber = 0; readNumber < readsPerChomosome; readNumber++) {
                out.addAlignment(createArtificialRead(header, "Read_" + readNumber, x - startingChromosome, readNumber, 100));
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
            for (int readNumber = 0; readNumber < readsPerChomosome; readNumber++) {
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
     * Create an artificial read based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param length         the length of the read
     * @return the artificial read
     */
    public static SAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceIndex(refIndex);
        record.setAlignmentStart(alignmentStart + 1);
        List<CigarElement> elements = new ArrayList<CigarElement>();
        elements.add(new CigarElement(length, CigarOperator.characterToEnum('M')));
        record.setCigar(new Cigar(elements));
        record.setProperPairFlag(false);
        return record;
    }


    /**
     * create an iterator containing the specified read piles
     * @param startingChr the chromosome (reference ID) to start from
     * @param endingChr the id to end with
     * @param readCount the number of reads per chromosome
     * @return StingSAMIterator representing the specified amount of fake data
     */
    public static StingSAMIterator unmappedReadIterator(int startingChr, int endingChr, int readCount ) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);
        Map<Integer, Integer> map = new HashMap<Integer, Integer>();
        for (int x = startingChr; x < endingChr; x++) {
            map.put(x,readCount);
        }
        return new ArtificialSAMIterator(startingChr,  endingChr,  readCount,header);
    }


}

/**
 * this fake iterator allows us to look at how specific piles of reads are handled
 */
class ArtificialSAMIterator implements StingSAMIterator {


    private int currentChromo = 0;
    private int currentRead = 0;
    private int readCount = 0;
    private boolean done = false;
    // the next record
    private SAMRecord next = null;
    SAMFileHeader header = null;

    // the passed in parameters
    final int sChr;
    final int eChr;
    final int rCount;

    /**
     * create the fake iterator, given the mapping of chromosomes and read counts
     * @param startingChr the starting chromosome
     * @param endingChr the ending chromosome
     * @param readCount the number of reads in each chromosome
     * @param header the associated header
     */
    ArtificialSAMIterator(int startingChr, int endingChr, int readCount, SAMFileHeader header) {
        sChr = startingChr;
        eChr = endingChr;
        rCount = readCount;
        this.header = header;
        this.currentChromo = startingChr;
    }

    public Reads getSourceInfo() {
        throw new UnsupportedOperationException("We don't support this");
    }

    public void close() {
        // done
        currentChromo = Integer.MAX_VALUE;
    }

    public boolean hasNext() {
        if (currentRead >= rCount) {
            currentChromo++;
            currentRead = 0;
        }
        if (currentChromo >= eChr) {
            return false;
        }
        int alignment = 0;
        if (currentChromo >= 0) {
            alignment = currentRead;
        } else {
            alignment = 0;
        }

        this.next = ArtificialSamUtils.createArtificialRead(this.header,String.valueOf(readCount), currentChromo, alignment, 50);
        currentRead++;
        return true;
    }

    public SAMRecord next() {
        return next;
    }

    public void remove() {
        throw new UnsupportedOperationException("You've tried to remove on a StingSAMIterator (unsupported), not to mention that this is a fake iterator.");
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
