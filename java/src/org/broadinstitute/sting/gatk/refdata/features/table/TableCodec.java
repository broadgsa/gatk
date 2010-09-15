package org.broadinstitute.sting.gatk.refdata.features.table;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.interval.IntervalUtils;

import java.io.IOException;
import java.util.*;

/**
 * implementation of a simple table (tab or comma delimited format) input files... more improvements to come
 */
public class TableCodec implements FeatureCodec  {
    private String delimiterRegex = "\\s+";
    private String headerDelimiter = "HEADER";
    private String commentDelimiter = "#";
    private ArrayList<String> header = new ArrayList<String>();

    @Override
    public Feature decodeLoc(String line) {
        return decode(line);
    }

    @Override
    public Feature decode(String line) {
        if (line.startsWith(headerDelimiter) || line.startsWith(commentDelimiter))
            return null;
        String[] split = line.split(delimiterRegex);
        if (split.length < 1)
            throw new IllegalArgumentException("TableCodec line = " + line + " doesn't appear to be a valid table format");


        return new TableFeature(GenomeLocParser.parseGenomeLoc(split[0]),Arrays.asList(split),header);
    }

    @Override
    public Class getFeatureType() {
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
