/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.codecs.table;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.gatk.utils.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

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
    protected final static String delimiterRegex = "\\s+";
    protected final static String headerDelimiter = "HEADER";
    protected final static String igvHeaderDelimiter = "track";
    protected final static String commentDelimiter = "#";
    // codec file extension
    protected final static String FILE_EXT = "tbl";

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
        return new TableFeature(genomeLocParser.parseGenomeLoc(split[0]),Arrays.asList(split), header);
    }

    /**
     * Can the file be decoded?
     * @param path path the file to test for parsability with this codec
     * @return true if the path has the correct file extension, false otherwise
     */
    @Override
    public boolean canDecode(final String path) { return path.endsWith("." + FILE_EXT); }

    @Override
    public Object readActualHeader(final LineIterator reader) {
        boolean isFirst = true;
        while (reader.hasNext()) {
            final String line = reader.peek(); // Peek to avoid reading non-header data
            if ( isFirst && ! line.startsWith(headerDelimiter) && ! line.startsWith(commentDelimiter)) {
                throw new UserException.MalformedFile("TableCodec file does not have a header");
            }
            isFirst &= line.startsWith(commentDelimiter);
            if (line.startsWith(headerDelimiter)) {
                reader.next(); // "Commit" the peek
                if (header.size() > 0) throw new IllegalStateException("Input table file seems to have two header lines.  The second is = " + line);
                final String spl[] = line.split(delimiterRegex);
                Collections.addAll(header, spl);
                return header;
            } else if (line.startsWith(commentDelimiter)) {
                reader.next(); // "Commit" the peek
            } else {
                break;
            }
        }
        return header;
    }
}
