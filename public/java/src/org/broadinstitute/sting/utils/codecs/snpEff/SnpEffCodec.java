/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.utils.codecs.snpEff;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.LineReader;
import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.EffectType;
import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.ChangeType;
import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.Zygosity;

import java.io.IOException;

public class SnpEffCodec implements FeatureCodec {

    public static final int EXPECTED_NUMBER_OF_FIELDS = 23;
    public static final String FIELD_DELIMITER_PATTERN = "\\t";
    public static final String EFFECT_FIELD_DELIMITER_PATTERN = "[,:]";
    public static final String HEADER_LINE_START = "# ";
    public static final String[] HEADER_FIELD_NAMES = { "Chromo",
                                                        "Position",
                                                        "Reference",
                                                        "Change",
                                                        "Change type",
                                                        "Homozygous",
                                                        "Quality",
                                                        "Coverage",
                                                        "Warnings",
                                                        "Gene_ID",
                                                        "Gene_name",
                                                        "Bio_type",
                                                        "Trancript_ID",   // yes, this is how it's spelled in the SnpEff output
                                                        "Exon_ID",
                                                        "Exon_Rank",
                                                        "Effect",
                                                        "old_AA/new_AA",
                                                        "Old_codon/New_codon",
                                                        "Codon_Num(CDS)",
                                                        "CDS_size",
                                                        "Codons around",
                                                        "AAs around",
                                                        "Custom_interval_ID"
                                                      };
    public static final int[] REQUIRED_FIELDS = { 0, 1, 15 };
    public static final String NON_CODING_GENE_FLAG = "WITHIN_NON_CODING_GENE";

    public Feature decodeLoc ( String line ) {
        return decode(line);
    }

    public Feature decode ( String line ) {
        String[] tokens = line.split(FIELD_DELIMITER_PATTERN, -1);

        if ( tokens.length != EXPECTED_NUMBER_OF_FIELDS ) {
            throw new TribbleException.InvalidDecodeLine("Line does not have the expected (" + EXPECTED_NUMBER_OF_FIELDS +
                                                         ") number of fields: found " + tokens.length + " fields.", line);
        }

        try {
            trimAllFields(tokens);
            checkForRequiredFields(tokens, line);

            String contig = tokens[0];
            long position = Long.parseLong(tokens[1]);

            String reference = tokens[2].isEmpty() ? null : tokens[2];
            String change = tokens[3].isEmpty() ? null : tokens[3];
            ChangeType changeType = tokens[4].isEmpty() ? null : ChangeType.valueOf(tokens[4]);
            Zygosity zygosity = tokens[5].isEmpty() ? null : Zygosity.valueOf(tokens[5]);
            Double quality = tokens[6].isEmpty() ? null : Double.parseDouble(tokens[6]);
            Long coverage = tokens[7].isEmpty() ? null : Long.parseLong(tokens[7]);
            String warnings = tokens[8].isEmpty() ? null : tokens[8];
            String geneID = tokens[9].isEmpty() ? null : tokens[9];
            String geneName = tokens[10].isEmpty() ? null : tokens[10];
            String bioType = tokens[11].isEmpty() ? null : tokens[11];
            String transcriptID = tokens[12].isEmpty() ? null : tokens[12];
            String exonID = tokens[13].isEmpty() ? null : tokens[13];
            Integer exonRank = tokens[14].isEmpty() ? null : Integer.parseInt(tokens[14]);

            boolean isNonCodingGene = isNonCodingGene(tokens[15]);
            int effectFieldTokenLimit = isNonCodingGene ? 3 : 2;
            String[] effectFieldTokens = tokens[15].split(EFFECT_FIELD_DELIMITER_PATTERN, effectFieldTokenLimit);
            EffectType effect = parseEffect(effectFieldTokens, isNonCodingGene);
            String effectExtraInformation = parseEffectExtraInformation(effectFieldTokens, isNonCodingGene);

            String oldAndNewAA = tokens[16].isEmpty() ? null : tokens[16];
            String oldAndNewCodon = tokens[17].isEmpty() ? null : tokens[17];
            Integer codonNum = tokens[18].isEmpty() ? null : Integer.parseInt(tokens[18]);
            Integer cdsSize = tokens[19].isEmpty() ? null : Integer.parseInt(tokens[19]);
            String codonsAround = tokens[20].isEmpty() ? null : tokens[20];
            String aasAround = tokens[21].isEmpty() ? null : tokens[21];
            String customIntervalID = tokens[22].isEmpty() ? null : tokens[22];

            return new SnpEffFeature(contig, position, reference, change, changeType, zygosity, quality, coverage,
                                     warnings, geneID, geneName, bioType, transcriptID, exonID, exonRank, isNonCodingGene,
                                     effect, effectExtraInformation, oldAndNewAA, oldAndNewCodon, codonNum, cdsSize,
                                     codonsAround, aasAround, customIntervalID);
        }
        catch ( NumberFormatException e ) {
            throw new TribbleException.InvalidDecodeLine("Error parsing a numeric field : " + e.getMessage(), line);
        }
        catch ( IllegalArgumentException e ) {
            throw new TribbleException.InvalidDecodeLine("Illegal value in field: " + e.getMessage(), line);
        }
    }

    private void trimAllFields ( String[] tokens ) {
        for ( int i = 0; i < tokens.length; i++ ) {
            tokens[i] = tokens[i].trim();
        }
    }

    private void checkForRequiredFields ( String[] tokens, String line ) {
        for ( int requiredFieldIndex : REQUIRED_FIELDS ) {
            if ( tokens[requiredFieldIndex].isEmpty() ) {
                throw new TribbleException.InvalidDecodeLine("Line is missing required field \"" +
                                                             HEADER_FIELD_NAMES[requiredFieldIndex] + "\"",
                                                             line);
            }
        }
    }

    private boolean isNonCodingGene ( String effectField ) {
        return effectField.startsWith(NON_CODING_GENE_FLAG);
    }

    private EffectType parseEffect ( String[] effectFieldTokens, boolean isNonCodingGene ) {
        String effectName = "";

        if ( effectFieldTokens.length > 1 && isNonCodingGene ) {
            effectName = effectFieldTokens[1].trim();
        }
        else {
            effectName = effectFieldTokens[0].trim();
        }

        return EffectType.valueOf(effectName);
    }

    private String parseEffectExtraInformation ( String[] effectFieldTokens, boolean isNonCodingGene ) {
        if ( (effectFieldTokens.length == 2 && ! isNonCodingGene) || effectFieldTokens.length == 3 ) {
            return effectFieldTokens[effectFieldTokens.length - 1].trim();
        }

        return null;
    }

    public Class getFeatureType() {
        return SnpEffFeature.class;
    }

    public Object readHeader ( LineReader reader ) {
        String headerLine = "";

        try {
            headerLine = reader.readLine();
        }
        catch ( IOException e ) {
            throw new TribbleException("Unable to read header line from input file.");
        }

        validateHeaderLine(headerLine);
        return headerLine;
    }

    private void validateHeaderLine ( String headerLine ) {
        if ( headerLine == null || ! headerLine.startsWith(HEADER_LINE_START) ) {
            throw new TribbleException.InvalidHeader("Header line does not start with " + HEADER_LINE_START);
        }

        String[] headerTokens = headerLine.substring(HEADER_LINE_START.length()).split(FIELD_DELIMITER_PATTERN);

        if ( headerTokens.length != EXPECTED_NUMBER_OF_FIELDS ) {
            throw new TribbleException.InvalidHeader("Header line does not contain headings for the expected number (" +
                                                     EXPECTED_NUMBER_OF_FIELDS + ") of columns.");
        }

        for ( int columnIndex = 0; columnIndex < headerTokens.length; columnIndex++ ) {
            if ( ! HEADER_FIELD_NAMES[columnIndex].equals(headerTokens[columnIndex]) ) {
                throw new TribbleException.InvalidHeader("Header field #" + columnIndex + ": Expected \"" +
                                                         HEADER_FIELD_NAMES[columnIndex] + "\" but found \"" +
                                                         headerTokens[columnIndex] + "\"");
            }
        }
    }
}
