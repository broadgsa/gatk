package org.broadinstitute.sting.gatk.refdata.features.sampileup;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.log4j.Logger;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.util.AsciiLineReader;
import org.broad.tribble.util.LineReader;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;

public class AnnotatorInputTableCodec implements FeatureCodec<AnnotatorInputTableFeature> {

    private static Logger logger = Logger.getLogger(AnnotatorInputTableCodec.class);

    public static final String DELIMITER = "\t";

    private ArrayList<String> header;

    /**
     * Parses the header.
     *
     * @param reader
     *
     * @return The # of header lines for this file.
     */
    public int readHeader(LineReader reader)
    {
        int[] lineCounter = new int[1];
        try {
            header = readHeader(reader, lineCounter);
        } catch(IOException e) {
            throw new IllegalArgumentException("Unable to read from file.", e);
        }
        return lineCounter[0];
    }


    /**
     * Parses the line into an AnnotatorInputTableFeature object.
     *
     * @param line
     */
    public AnnotatorInputTableFeature decode(String line) {
        final ArrayList<String> header = this.header; //optimization
        final ArrayList<String> values = Utils.split(line, DELIMITER, header.size());

        //if ( values.size() > header.size()) {
        //    throw new CodecLineParsingException(String.format("Encountered a line within " + file + " that has %d columns which is > the number of columns in the header which has %d columns.\nHeader: " + header + "\nLine: " + values, values.size(), header.size()));
        //}

        final AnnotatorInputTableFeature feature = new AnnotatorInputTableFeature(header);
        for ( int i = 0; i < header.size(); i++ ) {
            feature.putColumnValue(header.get(i), values.get(i));
        }

        final GenomeLoc loc = GenomeLocParser.parseGenomeLoc(values.get(0)); //GenomeLocParser.parseGenomeInterval(values.get(0)); - TODO switch to this

        //parse the location
        feature.setChr(loc.getContig());
        feature.setStart((int) loc.getStart());
        feature.setEnd((int) loc.getStop());

        return feature;
    }



    /**
     * Returns the header.
     * @param source
     * @return
     * @throws IOException
     */
    public static ArrayList<String> readHeader(final File source) throws IOException {
        FileInputStream is = new FileInputStream(source);
        try {
            return readHeader(new AsciiLineReader(is), null);
        } finally {
            is.close();
        }
    }


    /**
     * Returns the header, and also sets the 2nd arg to the number of lines in the header.
     * @param source
     * @param lineCounter An array of length 1 or null. If not null, array[0] will be set to the number of lines in the header.
     * @return The header fields.
     * @throws IOException
     */
    private static ArrayList<String> readHeader(final LineReader source, int[] lineCounter) throws IOException {

        ArrayList<String> header = null;
        int numLines = 0;

        //find the 1st line that's non-empty and not a comment
        String line = null;
        while( (line = source.readLine()) != null ) {
            numLines++;
            line = line.trim();
            if ( line.isEmpty() || line.startsWith("#") ) {
                continue;
            }

            //parse the header
            header = Utils.split(line, DELIMITER);
            break;
        }

        // check that we found the header
        if ( header == null ) {
            throw new IllegalArgumentException("No header in " + source + ". All lines are either comments or empty.");
        }

        if(lineCounter != null) {
            lineCounter[0] = numLines;
        }
        logger.info(String.format("Found header line containing %d columns:\n[%s]", header.size(), Utils.join("\t", header)));

        return header;
    }

}
