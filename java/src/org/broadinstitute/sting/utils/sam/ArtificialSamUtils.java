package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

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
    /**
     * create an artificial sam file
     * @param filename the filename to write to
     * @param numberOfChromosomes the number of chromosomes
     * @param startingChromosome where to start counting
     * @param chromosomeSize how large each chromosome is
     * @param readsPerChomosome how many reads to make in each chromosome.  They'll be aligned from position 1 to x (which is the number of reads)
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
     * @param filename the filename to write to
     * @param numberOfChromosomes the number of chromosomes
     * @param startingChromosome where to start counting
     * @param chromosomeSize how large each chromosome is
     * @param readsPerChomosome how many reads to make in each chromosome.  They'll be aligned from position 1 to x (which is the number of reads)
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
     * @param numberOfChromosomes the number of chromosomes to create
     * @param startingChromosome the starting number for the chromosome (most likely set to 1)
     * @param chromosomeSize the length of each chromosome
     * @return
     */
    public static SAMFileHeader createArtificialSamHeader(int numberOfChromosomes, int startingChromosome, int chromosomeSize) {
        SAMFileHeader header = new SAMFileHeader();
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        // make up some sequence records
        for (int x = startingChromosome; x < startingChromosome + numberOfChromosomes; x++) {
            SAMSequenceRecord rec = new SAMSequenceRecord("chr" + (x), chromosomeSize /* size */);
            rec.setSequenceLength(1000);
            dict.addSequence(rec);
        }
        header.setSequenceDictionary(dict);
        return header;
    }

    /**
     * Create an artificial read based on the parameters.  The cigar string will be *M, where * is the length of the read
     * @param header the SAM header to associate the read with
     * @param name the name of the read
     * @param refIndex the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param length the length of the read
     * @return the artificial read
     */
    public static SAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceIndex(refIndex);
        record.setAlignmentStart(alignmentStart + 1);
        List<CigarElement> elements = new ArrayList<CigarElement>();
        elements.add(new CigarElement( length, CigarOperator.characterToEnum('M')));
        record.setCigar(new Cigar(elements));
        record.setProperPairFlag(false);
        return record;
    }

}
