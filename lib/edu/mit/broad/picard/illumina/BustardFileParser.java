/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.illumina;

import edu.mit.broad.picard.util.PasteParser;
import edu.mit.broad.picard.util.FormatUtil;
import edu.mit.broad.picard.util.BasicTextFileParser;
import edu.mit.broad.picard.PicardException;

import java.io.File;
import java.io.FilenameFilter;
import java.io.Closeable;
import java.util.*;

/**
 * Class to parse the data in an Illumina Bustard directory and return an iterator over that data, in order
 * by tile.
 *
 * @author Kathleen Tibbetts
 */
public class BustardFileParser implements Iterator<BustardReadData>, Iterable<BustardReadData>, Closeable {

    private final File bustardDirectory;
    private final int lane;
    private final boolean pairedEnd;
    private PasteParser parser;
    private BustardReadData next = null;
    private final FormatUtil formatter = new FormatUtil();
    private boolean iterating = false;

    /**
     * Constructor
     *
     * @param bustardDirectory  directory where the Bustard files can be located
     * @param lane              the lane to parse
     * @param pairedEnd         whether this is a paired-end run
     */
    public BustardFileParser(File bustardDirectory, int lane, boolean pairedEnd) {
        this.bustardDirectory = bustardDirectory;
        this.lane = lane;
        this.pairedEnd = pairedEnd;
        initialize();
    }

    /**
     * Finds the relevant files in the bustardDirectory, sorts them, and puts them into the
     * <code>sortedFiles</code> iterator.  Does some basic sanity checking to ensure that some files
     * are found and that they are the expected multiple for paired-end or not.
     * 
     */
    private void initialize()
    {
        final String qseq1Regex = "s_" + lane + "_1_\\d{4}_qseq.txt(.gz)?";
        final String qseq2Regex = "s_" + lane + "_2_\\d{4}_qseq.txt(.gz)?";
        final String intensityRegex = "s_" + lane + "_\\d{4}_sig2.txt(.gz)?";

        File read1files[] = bustardDirectory.listFiles( new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.matches(qseq1Regex);
            }
        });

        File read2files[] = bustardDirectory.listFiles( new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.matches(qseq2Regex);
            }
        });

        File intensityFiles[] = bustardDirectory.listFiles( new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.matches(intensityRegex);
            }
        });

        // Some basic sanity checking on file counts
        if (read1files.length == 0 && read2files.length == 0 && intensityFiles.length == 0) {
            throw new PicardException("No Bustard files found in " +
                    bustardDirectory.getAbsolutePath() + " for lane " + lane);
        }
        if (pairedEnd) {
            if (read1files.length != read2files.length || read2files.length != intensityFiles.length) {
                throw new PicardException("Incorrect number of Bustard files found in " +
                        bustardDirectory.getAbsolutePath() + " for lane " + lane + ".  Found " +
                        read1files.length + " read 1 qseq files, " + read2files.length + " read 2 " +
                        "qseq files, and " + intensityFiles.length + " sig2 files.  There should be " +
                        "the same number of each type of file");
            }
        }
        else {
            if (read1files.length != intensityFiles.length) {
                throw new PicardException("Incorrect number of Bustard files found in " +
                        bustardDirectory.getAbsolutePath() + " for lane " + lane + ".  Found " +
                        read1files.length + " qseq files and " + intensityFiles.length + " sig2 files, " +
                        "which should be equal.");
            }
            if (read2files.length > 0) {
                throw new PicardException("Read 2 Bustard files found in " +
                        bustardDirectory.getAbsolutePath() + " for lane " + lane + ".  Lane " +
                        " was specified as a non-PE run, and so should not have any read 2 data.");
            }
        }

        // Sort each set of reads and create a text parser for it
        SortedSet<File> sortedRead1 = new TreeSet<File>(new BustardFilenameComparator());
        sortedRead1.addAll(Arrays.asList(read1files));
        read1files = sortedRead1.toArray(read1files);
        BasicTextFileParser read1Parser = new BasicTextFileParser(true, read1files);

        SortedSet<File> sortedIntensity = new TreeSet<File>(new BustardFilenameComparator());
        sortedIntensity.addAll(Arrays.asList(intensityFiles));
        intensityFiles = sortedIntensity.toArray(intensityFiles);
        BasicTextFileParser intensityParser = new BasicTextFileParser(true, intensityFiles);

        // And create a paste parser for all of them
        if (pairedEnd) {
            SortedSet<File> sortedRead2 = new TreeSet<File>(new BustardFilenameComparator());
            sortedRead2.addAll(Arrays.asList(read2files));
            read2files = sortedRead2.toArray(read2files);
            BasicTextFileParser read2Parser = new BasicTextFileParser(true, read2files);

            parser = new PasteParser(read1Parser, read2Parser, intensityParser);
        }
        else {
            parser = new PasteParser(read1Parser, intensityParser);
        }
    }

    /**
     * Parses the next line from the parser and constructs a BustardReadData object from it
     * The first 11 fields are the read1 data, the second 11 are the read2 data, and the remaining
     * values are the intensities data.  Note that the first four values in the intensity file
     * are not intensities but rather lane, tiles, x, and y for the given cluster.
     *
     * @param validate  whether to check that the expected number of intensity values are returned
     * @return  a fully populated BustardReadData object
     */
    private BustardReadData readNext(boolean validate) {
        if (!parser.hasNext()) {
                return null;
        }
        String data[][] = parser.next();
        String machine = data[0][0];
        int run = formatter.parseInt(data[0][1]);
        int lane = formatter.parseInt(data[0][2]);
        int tile = formatter.parseInt(data[0][3]);
        int x = formatter.parseInt(data[0][4]);
        int y = formatter.parseInt(data[0][5]);
        String firstSeq = data[0][8];
        String firstQual = data[0][9];
        boolean pf = formatter.parseInt(data[0][10]) == 1;
        String secondSeq = null;
        String secondQual = null;

        int intensityIndex = 1;
        if (pairedEnd) {
            secondSeq = data[1][8];
            secondQual = data[1][9];
            intensityIndex = 2;
        }

        int numIntensities = firstSeq.length() * (pairedEnd ? 2 : 1);

        // Sanity check since some of those files look a little weird
        if (validate) {
            int remaining = data[intensityIndex].length - 4;
            if ((remaining % 4 != 0) || (remaining/4) != numIntensities) {
                throw new PicardException("Unexpected number of intensity fields for " + machine + "/" + run +
                        "/" + lane + "/" + tile + ": " + remaining);
            }
        }

        double intensities[][] = new double[numIntensities][4];
        int intensityArrayIndex = 4;
        for (int i = 0; i < numIntensities; i++) {
            for (int j = 0; j < 4; j++) {
                intensities[i][j] = formatter.parseDouble(data[intensityIndex][intensityArrayIndex++]);
            }
        }

        return new BustardReadData(
                machine, run, lane, tile, firstSeq, firstQual, secondSeq, secondQual, pf, intensities, x, y);
        
    }

    /**
     * Returns an iterator over a set of elements of type BustardReadData.
     *
     * @return an iterator over a set of elements of type BustardReadData
     */
    public Iterator<BustardReadData> iterator() {
        if (iterating) {
            throw new IllegalStateException("iterator() method can only be called once, before the" +
                    "first call to hasNext()");
        }
        next = readNext(true);
        iterating = true;
        return this;
    }

    /**
     * Returns true if the iteration has more elements.
     *
     * @return  true if the iteration has more elements.  Otherwise returns false.
     */
    public boolean hasNext() {
        if (!iterating) {
            next = readNext(true);
            iterating = true;
        }
        return next != null;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return  the next element in the iteration
     * @throws java.util.NoSuchElementException
     */
    public BustardReadData next() {

        if (!hasNext()) {
            throw new NoSuchElementException("Iteration has no more elements.");
        }

        BustardReadData result = next;
        next = readNext(false);
        return result;
    }

    /**
     * Required method for Iterator API.
     *
     * @throws UnsupportedOperationException
     */
    public void remove() {
        throw new UnsupportedOperationException("Remove() not supported.");
    }

    /**
     * Closes the underlying PasteParser
     */
    public void close() {
        if (parser != null) {
            parser.close();
        }
    }

    public int getLane() { return this.lane; }
    public boolean isPairedEnd() { return this.pairedEnd; }
}
