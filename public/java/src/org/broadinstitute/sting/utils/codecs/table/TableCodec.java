/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.codecs.table;

import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Reads tab deliminated tabular text files
 *
 * <p>
 *     <ul>
 *     <li>Header: must begin with line HEADER or track (for IGV), followed by any number of column names,
 *     separated by whitespace.</li>
 *     <li>Comment lines starting with # are ignored</li>
 *     <li>Each non-header and non-comment line is split into parts by whitespace,
 *     and these parts are assigned as a map to their corresponding column name in the header.
 *     Note that the first element (corresponding to the HEADER column) must be a valid genome loc
 *     such as 1, 1:1 or 1:1-10, which is the position of the Table element on the genome.  TableCodec
 *     requires that there be one value for each column in the header, and no more, on all lines.</li>
 *     </ul>
 * </p>
 *
 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *     HEADER a b c
 *     1:1  1   2   3
 *     1:2  4   5   6
 *     1:3  7   8   9
 * </pre>
 *
 * @author Mark DePristo
 * @since 2009
 */
public class TableCodec extends AsciiFeatureCodec<TableFeature> implements ReferenceDependentFeatureCodec {
    final static protected String delimiterRegex = "\\s+";
    final static protected String headerDelimiter = "HEADER";
    final static protected String igvHeaderDelimiter = "track";
    final static protected String commentDelimiter = "#";

    protected ArrayList<String> header = new ArrayList<String>();

    /**
     * The parser to use when resolving genome-wide locations.
     */
    protected GenomeLocParser genomeLocParser;

    public TableCodec() {
        super(TableFeature.class);
    }

    /**
     * Set the parser to use when resolving genetic data.
     * @param genomeLocParser The supplied parser.
     */
    @Override
    public void setGenomeLocParser(GenomeLocParser genomeLocParser) {
        this.genomeLocParser =  genomeLocParser;
    }

    @Override
    public TableFeature decode(String line) {
        if (line.startsWith(headerDelimiter) || line.startsWith(commentDelimiter) || line.startsWith(igvHeaderDelimiter))
            return null;
        String[] split = line.split(delimiterRegex);
        if (split.length < 1)
            throw new IllegalArgumentException("TableCodec line = " + line + " doesn't appear to be a valid table format");
        return new TableFeature(genomeLocParser.parseGenomeLoc(split[0]),Arrays.asList(split),header);
    }

    @Override
    public Object readHeader(LineReader reader) {
        String line = "";
        try {
            boolean isFirst = true;
            while ((line = reader.readLine()) != null) {
                if ( isFirst && ! line.startsWith(headerDelimiter) && ! line.startsWith(commentDelimiter)) {
                    throw new UserException.MalformedFile("TableCodec file does not have a header");
                }
		isFirst &= line.startsWith(commentDelimiter);
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
