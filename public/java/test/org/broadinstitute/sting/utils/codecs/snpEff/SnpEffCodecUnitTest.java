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

import org.apache.commons.io.input.ReaderInputStream;
import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.EffectType;
import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.ChangeType;
import static org.broadinstitute.sting.utils.codecs.snpEff.SnpEffConstants.Zygosity;

import java.io.StringReader;

public class SnpEffCodecUnitTest {

    @Test
    public void testParseWellFormedSnpEffHeaderLine() {
        String wellFormedSnpEffHeaderLine = "# Chromo\tPosition\tReference\tChange\tChange type\t" +
                  "Homozygous\tQuality\tCoverage\tWarnings\tGene_ID\tGene_name\tBio_type\tTrancript_ID\tExon_ID\t" +
                  "Exon_Rank\tEffect\told_AA/new_AA\tOld_codon/New_codon\tCodon_Num(CDS)\tCDS_size\tCodons around\t" +
                  "AAs around\tCustom_interval_ID";

        SnpEffCodec codec = new SnpEffCodec();
        LineReader reader = new AsciiLineReader(new ReaderInputStream(new StringReader(wellFormedSnpEffHeaderLine)));
        String headerReturned = (String)codec.readHeader(reader);

        Assert.assertEquals(headerReturned, wellFormedSnpEffHeaderLine);
    }

    @Test(expectedExceptions = TribbleException.InvalidHeader.class)
    public void testParseWrongNumberOfFieldsSnpEffHeaderLine() {
        String wrongNumberOfFieldsSnpEffHeaderLine = "# Chromo\tPosition\tReference\tChange\tChange type\t" +
                  "Homozygous\tQuality\tCoverage\tWarnings\tGene_ID\tGene_name\tBio_type\tTrancript_ID\tExon_ID\t" +
                  "Exon_Rank\tEffect\told_AA/new_AA\tOld_codon/New_codon\tCodon_Num(CDS)\tCDS_size\tCodons around\t" +
                  "AAs around";

        SnpEffCodec codec = new SnpEffCodec();
        LineReader reader = new AsciiLineReader(new ReaderInputStream(new StringReader(wrongNumberOfFieldsSnpEffHeaderLine)));
        codec.readHeader(reader);
    }

    @Test(expectedExceptions = TribbleException.InvalidHeader.class)
    public void testParseMisnamedColumnSnpEffHeaderLine() {
        String misnamedColumnSnpEffHeaderLine = "# Chromo\tPosition\tRef\tChange\tChange type\t" +
                  "Homozygous\tQuality\tCoverage\tWarnings\tGene_ID\tGene_name\tBio_type\tTrancript_ID\tExon_ID\t" +
                  "Exon_Rank\tEffect\told_AA/new_AA\tOld_codon/New_codon\tCodon_Num(CDS)\tCDS_size\tCodons around\t" +
                  "AAs around\tCustom_interval_ID";

        SnpEffCodec codec = new SnpEffCodec();
        LineReader reader = new AsciiLineReader(new ReaderInputStream(new StringReader(misnamedColumnSnpEffHeaderLine)));
        codec.readHeader(reader);
    }

    @Test
    public void testParseSimpleEffectSnpEffLine() {
        String simpleEffectSnpEffLine = "1\t69428\tT\tG\tSNP\tHom\t6049.69\t61573\t\tENSG00000177693\t" +
                  "OR4F5\tmRNA\tENST00000326183\texon_1_69055_70108\t1\tNON_SYNONYMOUS_CODING\tF/C\tTTT/TGT\t113\t918\t\t\t";

        SnpEffFeature expectedFeature = new SnpEffFeature("1",
                                                          69428l,
                                                          "T",
                                                          "G",
                                                          ChangeType.SNP,
                                                          Zygosity.Hom,
                                                          6049.69,
                                                          61573l,
                                                          null,
                                                          "ENSG00000177693",
                                                          "OR4F5",
                                                          "mRNA",
                                                          "ENST00000326183",
                                                          "exon_1_69055_70108",
                                                          1,
                                                          false,
                                                          EffectType.NON_SYNONYMOUS_CODING,
                                                          null,
                                                          "F/C",
                                                          "TTT/TGT",
                                                          113,
                                                          918,
                                                          null,
                                                          null,
                                                          null
                                                         );

        SnpEffCodec codec = new SnpEffCodec();
        SnpEffFeature feature = (SnpEffFeature)codec.decode(simpleEffectSnpEffLine);

        Assert.assertEquals(feature, expectedFeature);
    }

    @Test
    public void testParseNonCodingRegionSnpEffLine() {
        String nonCodingRegionSnpEffLine = "1\t1337592\tG\tC\tSNP\tHom\t1935.52\t21885\t\tENSG00000250188\t" +
                  "RP4-758J18.5\tmRNA\tENST00000514958\texon_1_1337454_1338076\t2\tWITHIN_NON_CODING_GENE, NON_SYNONYMOUS_CODING\t" +
                  "L/V\tCTA/GTA\t272\t952\t\t\t";

        SnpEffFeature expectedFeature = new SnpEffFeature("1",
                                                          1337592l,
                                                          "G",
                                                          "C",
                                                          ChangeType.SNP,
                                                          Zygosity.Hom,
                                                          1935.52,
                                                          21885l,
                                                          null,
                                                          "ENSG00000250188",
                                                          "RP4-758J18.5",
                                                          "mRNA",
                                                          "ENST00000514958",
                                                          "exon_1_1337454_1338076",
                                                          2,
                                                          true,
                                                          EffectType.NON_SYNONYMOUS_CODING,
                                                          null,
                                                          "L/V",
                                                          "CTA/GTA",
                                                          272,
                                                          952,
                                                          null,
                                                          null,
                                                          null
                                                         );

        SnpEffCodec codec = new SnpEffCodec();
        SnpEffFeature feature = (SnpEffFeature)codec.decode(nonCodingRegionSnpEffLine);

        Assert.assertEquals(feature, expectedFeature);
    }

    @Test
    public void testParseExtraEffectInformationSnpEffLine() {
        String extraEffectInformationSnpEffLine = "1\t879537\tT\tC\tSNP\tHom\t341.58\t13733\t\tENSG00000187634\tSAMD11\t" +
                  "mRNA\tENST00000341065\t\t\tUTR_3_PRIME: 4 bases from transcript end\t\t\t\t\t\t\t";

        SnpEffFeature expectedFeature = new SnpEffFeature("1",
                                                          879537l,
                                                          "T",
                                                          "C",
                                                          ChangeType.SNP,
                                                          Zygosity.Hom,
                                                          341.58,
                                                          13733l,
                                                          null,
                                                          "ENSG00000187634",
                                                          "SAMD11",
                                                          "mRNA",
                                                          "ENST00000341065",
                                                          null,
                                                          null,
                                                          false,
                                                          EffectType.UTR_3_PRIME,
                                                          "4 bases from transcript end",
                                                          null,
                                                          null,
                                                          null,
                                                          null,
                                                          null,
                                                          null,
                                                          null
                                                         );

        SnpEffCodec codec = new SnpEffCodec();
        SnpEffFeature feature = (SnpEffFeature)codec.decode(extraEffectInformationSnpEffLine);

        Assert.assertEquals(feature, expectedFeature);
    }

    @Test
    public void testParseMultiEffectSnpEffLine() {
        String multiEffectSnpEffLine = "1\t901901\tC\tT\tSNP\tHom\t162.91\t4646\t\tENSG00000187583\tPLEKHN1\tmRNA\t" +
                  "ENST00000379410\texon_1_901877_901994\t1\tSTART_GAINED: ATG, UTR_5_PRIME: 11 bases from TSS\t\t\t\t\t\t\t";

        SnpEffFeature expectedFeature = new SnpEffFeature("1",
                                                          901901l,
                                                          "C",
                                                          "T",
                                                          ChangeType.SNP,
                                                          Zygosity.Hom,
                                                          162.91,
                                                          4646l,
                                                          null,
                                                          "ENSG00000187583",
                                                          "PLEKHN1",
                                                          "mRNA",
                                                          "ENST00000379410",
                                                          "exon_1_901877_901994",
                                                          1,
                                                          false,
                                                          EffectType.START_GAINED,
                                                          "ATG, UTR_5_PRIME: 11 bases from TSS",
                                                          null,
                                                          null,
                                                          null,
                                                          null,
                                                          null,
                                                          null,
                                                          null
                                                         );

        SnpEffCodec codec = new SnpEffCodec();
        SnpEffFeature feature = (SnpEffFeature)codec.decode(multiEffectSnpEffLine);

        Assert.assertEquals(feature, expectedFeature);
    }

    @Test(expectedExceptions = TribbleException.InvalidDecodeLine.class)
    public void testParseWrongNumberOfFieldsSnpEffLine() {
        String wrongNumberOfFieldsSnpEffLine = "1\t69428\tT\tG\tSNP\tHom\t6049.69\t61573\t\tENSG00000177693\t" +
                  "OR4F5\tmRNA\tENST00000326183\texon_1_69055_70108\t1\tNON_SYNONYMOUS_CODING\tF/C\tTTT/TGT\t113\t918\t\t";

        SnpEffCodec codec = new SnpEffCodec();
        SnpEffFeature feature = (SnpEffFeature)codec.decode(wrongNumberOfFieldsSnpEffLine);
    }

    @Test(expectedExceptions = TribbleException.InvalidDecodeLine.class)
    public void testParseBlankEffectFieldSnpEffLine() {
        String blankEffectFieldSnpEffLine = "1\t69428\tT\tG\tSNP\tHom\t6049.69\t61573\t\tENSG00000177693\t" +
                  "OR4F5\tmRNA\tENST00000326183\texon_1_69055_70108\t1\t\tF/C\tTTT/TGT\t113\t918\t\t\t";

        SnpEffCodec codec = new SnpEffCodec();
        SnpEffFeature feature = (SnpEffFeature)codec.decode(blankEffectFieldSnpEffLine);
    }

    @Test(expectedExceptions = TribbleException.InvalidDecodeLine.class)
    public void testParseInvalidNumericFieldSnpEffLine() {
        String invalidNumericFieldSnpEffLine = "1\t69428\tT\tG\tSNP\tHom\t6049.69\t61573\t\tENSG00000177693\t" +
                  "OR4F5\tmRNA\tENST00000326183\texon_1_69055_70108\t1\tNON_SYNONYMOUS_CODING\tF/C\tTTT/TGT\t113\tfoo\t\t\t";;

        SnpEffCodec codec = new SnpEffCodec();
        SnpEffFeature feature = (SnpEffFeature)codec.decode(invalidNumericFieldSnpEffLine);
    }
}
