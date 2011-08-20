package org.broadinstitute.sting.utils.codecs.beagle;
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


import org.broad.tribble.Feature;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

/**
 * TODO GUILLERMO DEL ANGEL
 *
 * <p>
 * Codec Description
 * </p>
 *
 * <p>
 * See also: @see <a href="http://vcftools.sourceforge.net/specs.html">VCF specification</a><br>
 * </p>

 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *     line 1
 *     line 2
 *     line 3
 * </pre>
 *
 * @author Mark DePristo
 * @since 2010
 */
public class BeagleCodec implements ReferenceDependentFeatureCodec<BeagleFeature> {
    private String[] header;
    public enum BeagleReaderType {PROBLIKELIHOOD, GENOTYPES, R2};
    private BeagleReaderType readerType;
    private int valuesPerSample;
    private int initialSampleIndex;
    private int markerPosition;
    private ArrayList<String> sampleNames;
    private int expectedTokensPerLine;

    private static final String delimiterRegex = "\\s+";

    /**
     * The parser to use when resolving genome-wide locations.
     */
    private GenomeLocParser genomeLocParser;

    /**
     * Set the parser to use when resolving genetic data.
     * @param genomeLocParser The supplied parser.
     */
    public void setGenomeLocParser(GenomeLocParser genomeLocParser) {
        this.genomeLocParser =  genomeLocParser;
    }    

    public Feature decodeLoc(String line) {
        return decode(line);
    }
    
    public static String[] readHeader(final File source) throws IOException {
        FileInputStream is = new FileInputStream(source);
        try {
            return readHeader(new AsciiLineReader(is), null);
        } finally {
            is.close();
        }
    }

    public Object readHeader(LineReader reader)
    {
        int[] lineCounter = new int[1];
        try {
            header = readHeader(reader, lineCounter);

            Boolean getSamples = true;
            markerPosition = 0; //default value for all readers

            if (header[0].matches("I")) {
                // Phased genotype Beagle files start with "I"
                readerType = BeagleReaderType.GENOTYPES;
                valuesPerSample = 2;
                initialSampleIndex = 2;
                markerPosition = 1;
            }
            else if (header[0].matches("marker")) {
                readerType = BeagleReaderType.PROBLIKELIHOOD;
                valuesPerSample = 3;
                initialSampleIndex = 3;
            }
            else {
                readerType = BeagleReaderType.R2;
                getSamples = false;
                // signal we don't have a header
                lineCounter[0] = 0;
                // not needed, but for consistency:
                valuesPerSample = 0;
                initialSampleIndex = 0;
            }

            sampleNames = new ArrayList<String>();

            if (getSamples) {
                for (int k = initialSampleIndex; k < header.length; k += valuesPerSample)
                    sampleNames.add(header[k]);

                expectedTokensPerLine = sampleNames.size()*valuesPerSample+initialSampleIndex;

            } else {
                expectedTokensPerLine = 2;
            }


        } catch(IOException e) {
            throw new IllegalArgumentException("Unable to read from file.", e);
        }
        return header;
    }

    private static String[] readHeader(final LineReader source, int[] lineCounter) throws IOException {

        String[] header = null;
        int numLines = 0;

        //find the 1st line that's non-empty and not a comment
        String line;
        while( (line = source.readLine()) != null ) {
            numLines++;
            if ( line.trim().isEmpty() ) {
                continue;
            }

            //parse the header
            header = line.split(delimiterRegex);
            break;
        }

        // check that we found the header
        if ( header == null ) {
            throw new IllegalArgumentException("No header in " + source);
        }

        if(lineCounter != null) {
            lineCounter[0] = numLines;
        }

        return header;
    }

    private static Pattern MARKER_PATTERN = Pattern.compile("(.+):([0-9]+)");

    @Override
    public Class<BeagleFeature> getFeatureType() {
        return BeagleFeature.class;
    }

    public BeagleFeature decode(String line) {
        String[] tokens;

        // split the line
        tokens = line.split(delimiterRegex);
        if (tokens.length != expectedTokensPerLine)
            throw new CodecLineParsingException("Incorrect number of fields in Beagle input on line "+line);



        BeagleFeature bglFeature = new BeagleFeature();

        final GenomeLoc loc = genomeLocParser.parseGenomeLoc(tokens[markerPosition]); //GenomeLocParser.parseGenomeLoc(values.get(0)); - TODO switch to this

        //parse the location: common to all readers
        bglFeature.setChr(loc.getContig());
        bglFeature.setStart((int) loc.getStart());
        bglFeature.setEnd((int) loc.getStop());

        // Parse R2 if needed
        if (readerType == BeagleReaderType.R2) {
            bglFeature.setR2value(Double.valueOf(tokens[1]));
        }
        else if (readerType == BeagleReaderType.GENOTYPES) {
            // read phased Genotype pairs
            HashMap<String, ArrayList<String>> sampleGenotypes = new HashMap<String, ArrayList<String>>();

            for ( int i = 2; i < tokens.length; i+=2 ) {
                String sampleName = sampleNames.get(i/2-1);
                if ( ! sampleGenotypes.containsKey(sampleName) ) {
                    sampleGenotypes.put(sampleName, new ArrayList<String>());
                }

                sampleGenotypes.get(sampleName).add(tokens[i]);
                sampleGenotypes.get(sampleName).add(tokens[i+1]);
            }

            bglFeature.setGenotypes(sampleGenotypes);
        }
        else {
            // read probabilities/likelihood trios and alleles
            bglFeature.setAlleleA(tokens[1], true);
            bglFeature.setAlleleB(tokens[2], false);
            HashMap<String, ArrayList<String>> sampleProbLikelihoods = new HashMap<String, ArrayList<String>>();

            for ( int i = 3; i < tokens.length; i+=3 ) {
                String sampleName = sampleNames.get(i/3-1);
                if ( ! sampleProbLikelihoods.containsKey(sampleName) ) {
                    sampleProbLikelihoods.put(sampleName, new ArrayList<String>());
                }

                sampleProbLikelihoods.get(sampleName).add(tokens[i]);
                sampleProbLikelihoods.get(sampleName).add(tokens[i+1]);
                sampleProbLikelihoods.get(sampleName).add(tokens[i+2]);
            }
            bglFeature.setProbLikelihoods(sampleProbLikelihoods);
        }

        return bglFeature;
    }

}
