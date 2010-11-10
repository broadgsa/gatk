package org.broadinstitute.sting.utils.bed;

import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Oct 5, 2009
 * Time: 5:46:45 PM
 */
public class BedParser {
    // the GATk operates as a one based location, bed files are 0 based
    static final int TO_ONE_BASED_ADDITION = 1;

    // the buffered reader input
    private final BufferedReader mIn;

    private GenomeLocParser genomeLocParser;

    // our array of locations
    private List<GenomeLoc> mLocations;

    /**
     * parse a bed file, given it's location
     *
     * @param fl
     */
    public BedParser(GenomeLocParser genomeLocParser,File fl) {
        this.genomeLocParser = genomeLocParser;
        try {
            mIn = new BufferedReader(new FileReader(fl));
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(fl, e);
        }
        mLocations = parseLocations();
    }

    /**
     * parse a bed file, given an input reader
     *
     * @param fl the bed file
     */
    public BedParser(BufferedReader fl) {
        mIn = fl;
        mLocations = parseLocations();
    }

    /**
     * parse out the locations
     *
     * @return a list of GenomeLocs, sorted and merged
     */
    private List<GenomeLoc> parseLocations() {
        String line = null;
        List<GenomeLoc> locArray = new ArrayList<GenomeLoc>();
        try {
            while ((line = mIn.readLine()) != null) {
                locArray.add(parseLocation(genomeLocParser,line));
            }
        } catch (IOException e) {
            throw new UserException.MalformedFile("Unable to parse line in BED file.");
        }
        return locArray;
    }

    /**
     * parse a single location
     *
     * @param line the line, as a string
     * @return a parsed genome loc
     */
    public static GenomeLoc parseLocation(GenomeLocParser genomeLocParser,String line) {
        String contig;
        int start;
        int stop;
        try {
            String parts[] = line.split("\\s+");
            contig = parts[0];
            start = Integer.valueOf(parts[1]) + TO_ONE_BASED_ADDITION;
            stop = Integer.valueOf(parts[2]); // the ending point is an open interval
        } catch (Exception e) {
            throw new UserException.MalformedFile("Unable to process bed file line = " + line, e);
        }

        // we currently drop the rest of the bed record, which can contain names, scores, etc
        return genomeLocParser.createGenomeLoc(contig, start, stop);

    }

    /**
     * return the sorted, and merged (for overlapping regions)
     *
     * @return an arraylist
     */
    public List<GenomeLoc> getLocations() {
        return mLocations;
    }
}
