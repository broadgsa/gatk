/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.codecs.samread;

import net.sf.samtools.Cigar;
import net.sf.samtools.TextCigarCodec;
import net.sf.samtools.util.StringUtil;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.ParsingUtils;

/**
 * Decodes a simple SAM text string.
 *
 * <p>
 * Reads in the SAM text version of a BAM file as a ROD.  For testing only
 * </p>
 *
 * <p>
 * See also: @see <a href="http://samtools.sourceforge.net">SAMTools</a> for format specification
 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *     SL-XBC:1:10:628:923#0	16	Escherichia_coli_K12	1	37	76M	=	1	0	AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGA	B@>87<;A@?@957:>>@AA@B>@A9AB@B>@A@@@@@A;=AAB@BBBBBCBBBB@>A>:ABB@BAABCB=CA@CB
 * </pre>
 *
 * @author Matt Hanna
 * @since 2009
 */
public class SAMReadCodec implements FeatureCodec<SAMReadFeature> {
    /* SL-XBC:1:10:628:923#0	16	Escherichia_coli_K12	1	37	76M	=	1	0	AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGA	B@>87<;A@?@957:>>@AA@B>@A9AB@B>@A@@@@@A;=AAB@BBBBBCBBBB@>A>:ABB@BAABCB=CA@CB */

    // the number of tokens we expect to parse from a read line
    private static final int expectedTokenCount = 11;

    /**
     * Return the # of header lines for this file.
     *
     * @param reader the line reader
     * @return 0 in this case, we assume no header lines.  The reads file may have a
     *         header line beginning with '@', but we can ignore that in the decode function.
     */
    public Object readHeader(LineReader reader) {
        // we don't require a header line, but it may exist.  We'll deal with that above.
        return null;
    }

    @Override
    public Class<SAMReadFeature> getFeatureType() {
        return SAMReadFeature.class;
    }

    public Feature decodeLoc(String line) {
        return decode(line);
    }

    /**
     * Decode a single line in a SAM text file.
     * @param line line to decode.
     * @return A SAMReadFeature modeling that line.
     */
    public SAMReadFeature decode(String line) {
        // we may be asked to process a header line; ignore it
        if (line.startsWith("@")) return null;        

        String[] tokens = new String[expectedTokenCount];

        // split the line
        int count = ParsingUtils.splitWhitespace(line,tokens);

        // check to see if we've parsed the string into the right number of tokens (expectedTokenCount)
        if (count != expectedTokenCount)
            throw new CodecLineParsingException("the SAM read line didn't have the expected number of tokens " +
                                                "(expected = " + expectedTokenCount + ", saw = " + count + " on " +
                                                "line = " + line + ")");

        final String readName = tokens[0];
        final int flags = Integer.parseInt(tokens[1]);
        final String contigName = tokens[2];
        final int alignmentStart = Integer.parseInt(tokens[3]);
        final int mapQ = Integer.parseInt(tokens[4]);
        final String cigarString = tokens[5];
        final String mateContigName = tokens[6];
        final int mateAlignmentStart = Integer.parseInt(tokens[7]);
        final int inferredInsertSize = Integer.parseInt(tokens[8]);
        final byte[] bases = StringUtil.stringToBytes(tokens[9]);
        final byte[] qualities = StringUtil.stringToBytes(tokens[10]);

        // Infer the alignment end.
        Cigar cigar = TextCigarCodec.getSingleton().decode(cigarString);
        int alignmentEnd = alignmentStart + cigar.getReferenceLength() - 1;

        // Remove printable character conversion from the qualities.
        for(byte quality: qualities) quality -= 33;

        return new SAMReadFeature(readName,
                                  flags,
                                  contigName,
                                  alignmentStart,
                                  alignmentEnd,
                                  mapQ,
                                  cigarString,
                                  mateContigName,
                                  mateAlignmentStart,
                                  inferredInsertSize,
                                  bases,
                                  qualities);
    }


}
