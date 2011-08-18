package org.broadinstitute.sting.gatk.refdata.features.table;

import org.broad.tribble.Feature;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * implementation of a simple table (tab or comma delimited format) input files
 */
public class TableCodec implements ReferenceDependentFeatureCodec {
    protected String delimiterRegex = "\\s+";
    protected String headerDelimiter = "HEADER";
    protected String igvHeaderDelimiter = "track";
    protected String commentDelimiter = "#";
    protected ArrayList<String> header = new ArrayList<String>();

    /**
     * The parser to use when resolving genome-wide locations.
     */
    protected GenomeLocParser genomeLocParser;

    /**
     * Set the parser to use when resolving genetic data.
     * @param genomeLocParser The supplied parser.
     */
    @Override
    public void setGenomeLocParser(GenomeLocParser genomeLocParser) {
        this.genomeLocParser =  genomeLocParser;
    }


    @Override
    public Feature decodeLoc(String line) {
        return decode(line);
    }

    @Override
    public Feature decode(String line) {
        if (line.startsWith(headerDelimiter) || line.startsWith(commentDelimiter) || line.startsWith(igvHeaderDelimiter))
            return null;
        String[] split = line.split(delimiterRegex);
        if (split.length < 1)
            throw new IllegalArgumentException("TableCodec line = " + line + " doesn't appear to be a valid table format");
        return new TableFeature(genomeLocParser.parseGenomeLoc(split[0]),Arrays.asList(split),header);
    }

    @Override
    public Class<TableFeature> getFeatureType() {
        return TableFeature.class;
    }

    @Override
    public Object readHeader(LineReader reader) {
        String line = "";
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(headerDelimiter)) {
                    if (header.size() > 0) throw new IllegalStateException("Input table file seems to have two header lines.  The second is = " + line);
                    String spl[] = line.split(delimiterRegex);
                    for (String s : spl) header.add(s);
                    return header;
                } else if (!line.startsWith(commentDelimiter)) {
                    break;
                }
            }
        } catch (IOException e) {
            throw new UserException.MalformedFile("unable to parse header from TableCodec file",e);
        }
        return header;
    }
}
