/*
 * Copyright (c) 2010 The Broad Institute
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

package net.sf.picard.reference;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import static net.sf.picard.reference.FastaSequenceIndexBuilder.Status.*;

import java.io.*;

import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Builds FastaSequenceIndex from fasta file.
 * Produces fai file with same output as samtools faidx
 */
public class FastaSequenceIndexBuilder {
    final public File fastaFile;
    final boolean printProgress;

    // keep track of location in file
    long bytesRead, endOfLastLine, lastTimestamp, fileLength;  // initialized to -1 to keep 0-indexed position in file;

    // vars that store information about the contig that is currently being read
    String contig;
    long location, size, bytesPerLine, basesPerLine, basesThisLine;
    int thisSequenceIndex = 0;

    // vars that keep loop state
    byte lastByte = 0, currentByte = 0, nextByte = 0;
    public enum Status { NONE, CONTIG, FIRST_SEQ_LINE, SEQ_LINE, COMMENT }
    Status status = Status.NONE; // keeps state of what is currently being read. better to use int instead of enum?

    public FastaSequenceIndexBuilder(File fastaFile, boolean printProgress) {
        this.fastaFile = fastaFile;
        fileLength = fastaFile.length();
        this.printProgress = printProgress;
    }

    /**
     * Creates fasta sequence index from fasta file
     * @return FastaSequenceIndex that is read from file
     */
    public FastaSequenceIndex createIndex() {    // should this be static?
        bytesRead = -1;
        endOfLastLine = -1;
        contig = "";
        location = 0;
        size = 0;
        bytesPerLine = 0;
        basesPerLine = 0;
        basesThisLine = 0;
        lastTimestamp = System.currentTimeMillis();
        FastaSequenceIndex sequenceIndex = new FastaSequenceIndex();

        // initialize input stream
        DataInputStream in;
        try {
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(fastaFile)));
        }
        catch (Exception e) {
            throw new UserException.CouldNotReadInputFile(fastaFile, "Could not read fasta file", e);
        }

        /*
        * iterate through each character in file one at a time, but must account for variance in line terminators
        * strategy is to check if current character is a line terminator (\n, \r), then check next character
        * only treat as end of line if next character is NOT a line terminator
        */
        try {
            // intialize iterators
            nextByte = in.readByte();
            currentByte = '\n';
            while(currentByte != -1) {

                bytesRead ++; // update position in file
                lastByte = currentByte;
                currentByte = nextByte;
                try {
                    nextByte = in.readByte();
                }                                          
                catch (EOFException e) {
                    nextByte = -1;
                }

                switch(status) {

                    // if not currently reading anything
                    // only thing that can trigger action is '>' (start of contig) or ';' (comment)
                    case NONE:
                        if (currentByte == '>')
                            status = CONTIG;
                        else if (currentByte == ';')
                            status = COMMENT;
                        break;

                    // if reading a comment, just ignore everything until end of line
                    case COMMENT:
                        if (isEol(currentByte)) {
                            if (!isEol(nextByte))
                                status = Status.NONE;
                        }
                        break;

                    // if reading a contig, add char to contig string
                    // contig string can be anything, including spaces
                    case CONTIG:
                        if (isEol(currentByte)) {
                            if (!isEol(nextByte)) {
                                status = Status.FIRST_SEQ_LINE;
                                location = bytesRead + 1;
                            }
                        }
                        else
                            contig += (char) currentByte;
                        break;

                    // record bases and bytes of first sequence line, to validate against later lines
                    case FIRST_SEQ_LINE:

                        if (isEol(currentByte)) {

                            // record bases per line if last character was a base
                            if (!isEol(lastByte)) {
                                basesPerLine = bytesRead - location;
                                basesThisLine = basesPerLine;
                                size += basesPerLine;
                            }

                            // next character is start of next line, now know bytes per line
                            if (!isEol(nextByte)) {   // figure out what to do if there is only one data line
                                bytesPerLine = bytesRead - location + 1;
                                status = Status.SEQ_LINE;
                                endOfLastLine = bytesRead;

                                // if next char is ';' or '>', then there is only one contig =>
                                if (nextByte == ';' || nextByte == '>')
                                    finishReadingContig(sequenceIndex);
                            }
                        }

                        // validate base character
                        else {
                            if (!isValidBase(currentByte))
                                throw new UserException.MalformedFile(fastaFile, String.format("An invalid base was found in the contig: %s", contig));
                        }
                        break;


                    case SEQ_LINE:

                        if (isEol(currentByte)) {

                            // record bases per line if last character was a base
                            if (!isEol(lastByte)) {
                                basesThisLine = bytesRead - endOfLastLine - 1;
                                size += basesThisLine;
                            }

                            // reached end of line - check if end of contig
                            if (!isEol(nextByte)) {

                                // if comment or new contig, definitely end of sequence
                                if (nextByte == ';' || nextByte == '>')
                                    finishReadingContig(sequenceIndex);

                                    // if this line has different # of bases OR same # of bases and different # of bytes:
                                    // error if next char is a valid base, end of contig otherwise
                                else if (basesThisLine != basesPerLine || bytesPerLine != bytesRead - endOfLastLine) {
                                    if (isValidBase(nextByte) && nextByte != -1) {
                                        throw new UserException.MalformedFile(fastaFile, String.format("An invalid line was found in the contig: %s", contig));
                                    }
                                    else
                                        finishReadingContig(sequenceIndex);
                                }
                                endOfLastLine = bytesRead;
                            }
                        }

                        // validate base character
                        else {
                            if (!isValidBase(currentByte))
                                throw new UserException.MalformedFile(fastaFile, String.format("An invalid base was found in the contig: %s", contig));
                        }
                        break;
                }
            }
            return sequenceIndex;
        }
        catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(fastaFile, "Could not read fasta file", e);
        }
        catch (Exception e) {
            throw new ReviewedStingException(e.getMessage(), e);
        }
    }

    /**
     * Checks if character is an end of line character, \n or \r
     * @param currentByte Character to check, as a byte
     * @return True if current character is \n or \r' false otherwise
     */
    private boolean isEol(byte currentByte) {
        return (currentByte == '\n' || currentByte == '\r');
    }


    /**
     * checks that character is valid base
     * only checks that the base isn't whitespace, like picard does
     * could add check that character is A/C/G/T/U if wanted
     * @param currentByte Character to check, as a byte
     * @return True if character is not whitespace, false otherwise
     */
    private boolean isValidBase(byte currentByte) {
        return (!Character.isWhitespace(currentByte) && currentByte != ';' && currentByte != '>');
    }

    /*
     * When reader reaches the end of a contig
     * Reset iterators and add contig to sequence index
     */
    private void finishReadingContig(FastaSequenceIndex sequenceIndex) {
        sequenceIndex.add(new FastaSequenceIndexEntry(contig, location, size, (int) basesPerLine, (int) bytesPerLine, thisSequenceIndex++));
        status = Status.NONE;
        contig = "";
        size = 0;

        if (System.currentTimeMillis() - lastTimestamp > 10000) {
            int percentProgress = (int) (100*bytesRead/fileLength);
            if (printProgress)
                System.out.println(String.format("PROGRESS UPDATE: file is %d percent complete", percentProgress));
            lastTimestamp = System.currentTimeMillis();
        }
    }

    /**
     * Stores FastaSequenceIndex as a .fasta.fai file on local machine
     * Although method is public it cannot be called on any old FastaSequenceIndex - must be created by a FastaSequenceIndexBuilder
     * @param sequenceIndex sequenceIndex to be saved
     * @param faiFile file where we should store index
     */
    public static void saveAsFaiFile(FastaSequenceIndex sequenceIndex, File faiFile) {

        BufferedWriter out;
        try {
            out = new BufferedWriter(new FileWriter(faiFile));
        }
        catch (Exception e) {
            throw new UserException.CouldNotCreateOutputFile(faiFile, e);
        }

        try {
            for(FastaSequenceIndexEntry entry: sequenceIndex) {
                out.write(toIndexFileLine(entry));
                out.newLine();
            }
            out.close();
        }
        catch (Exception e) {
            throw new UserException.CouldNotCreateOutputFile(faiFile, e);
        }
    }

    /**
     * Print string in format of fai file line
     * @return Contig as one line in a fai file
     */
    private static String toIndexFileLine(FastaSequenceIndexEntry entry) {
        return String.format("%s\t%d\t%d\t%d\t%d", entry.getContig(), entry.getSize(), entry.getLocation(), entry.getBasesPerLine(), entry.getBytesPerLine());
    }

}
