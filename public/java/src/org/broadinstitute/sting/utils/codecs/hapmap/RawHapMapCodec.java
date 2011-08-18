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

package org.broadinstitute.sting.utils.codecs.hapmap;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;
import java.util.Arrays;

/**
 * a codec for the file types produced by the HapMap consortium, available on their website:
 * http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/
 *
 * The format includes eleven standard fields, plus genotypes for each of the samples included
 * in the file
 * 
 */
public class RawHapMapCodec implements FeatureCodec {
    // the minimum number of features in the HapMap file line
    private static final int minimumFeatureCount = 11;

    private String headerLine;
    /**
     * decode the location only
     * @param line the input line to decode
     * @return a HapMapFeature
     */
    public Feature decodeLoc(String line) {
        return decode(line);
    }

    /**
     * decode the hapmap record
     * @param line the input line to decode
     * @return a HapMapFeature, with the given fields 
     */
    public Feature decode(String line) {
        String[] array = line.split("\\s+");

        // make sure the split was successful - that we got an appropriate number of fields
        if (array.length < minimumFeatureCount)
            throw new IllegalArgumentException("Unable to parse line " + line + ", the length of split features is less than the minimum of " + minimumFeatureCount);

        // create a new feature given the array
        return new RawHapMapFeature(array[0],
                array[1].split("/"),
                array[2],
                Long.valueOf(array[3]),
                Strand.toStrand(array[4]),
                array[5],
                array[6],
                array[7],
                array[8],
                array[9],
                array[10],
                Arrays.copyOfRange(array,11,array.length),
                headerLine);
    }

    public Class<RawHapMapFeature> getFeatureType() {
        return RawHapMapFeature.class;
    }

    public Object readHeader(LineReader reader) {
        try {
            headerLine = reader.readLine();
        } catch (IOException e) {
            throw new IllegalArgumentException("Unable to read a line from the line reader");
        }
        return headerLine;
    }
}
