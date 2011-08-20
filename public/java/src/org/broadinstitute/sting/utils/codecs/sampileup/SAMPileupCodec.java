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

package org.broadinstitute.sting.utils.codecs.sampileup;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.ParsingUtils;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.broadinstitute.sting.utils.codecs.sampileup.SAMPileupFeature.VariantType;

/**
 * Decoder for SAM pileup data.  For GATK validation purposes only
 *
 * <p>
 *     Pileup format is first used by Tony Cox and Zemin Ning at the Sanger Institute.
 *     It desribes the base-pair information at each chromosomal position. This format
 *     facilitates SNP/indel calling and brief alignment viewing by eyes.
 * </p>
 * <p>
 *     Each line consists of chromosome, 1-based coordinate, reference base, the
 *     number of reads covering the site, read bases and base qualities. At the
 *     read base column, a dot stands for a match to the reference base on the
 *     forward strand, a comma for a match on the reverse strand, `ACGTN' for a mismatch
 *     on the forward strand and `acgtn' for a mismatch on the reverse strand.
 *     A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between
 *     this reference position and the next reference position. The length of the
 *     insertion is given by the integer in the pattern, followed by the inserted sequence.
 * </p>
 *
 * <p>
 *     <br>See also: @see <a href="http://samtools.sourceforge.net/">SAMTools project</a></br>
 *     <br>See also: @see <a href="http://samtools.sourceforge.net/pileup.shtml">Pileup format</a></br>
 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *     seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
 *     seq1 273 T 23  ,.....,,.,.,...,,,.,..A <<<;<<<<<<<<<3<=<<<;<<+
 *     seq1 274 T 23  ,.$....,,.,.,...,,,.,...    7<7;<;<<<<<<<<<=<;<;<<6
 *     seq1 275 A 23  ,$....,,.,.,...,,,.,...^l.  <+;9*<<<<<<<<<=<<:;<<<<
 *     seq1 276 G 22  ...T,,.,.,...,,,.,....  33;+<<7=7<<7<&<<1;<<6<
 *     seq1 277 T 22  ....,,.,.,.C.,,,.,..G.  +7<;<<<<<<<&<=<<:;<<&<
 *     seq1 278 G 23  ....,,.,.,...,,,.,....^k.   %38*<<;<7<<7<=<<<;<<<<<
 *     seq1 279 C 23  A..T,,.,.,...,,,.,..... ;75&<<<<<<<<<=<<<9<<:<<
 * </pre>
 *
 * @author Matt Hanna
 * @since 2009
 */
public class SAMPileupCodec implements FeatureCodec<SAMPileupFeature> {
    // the number of tokens we expect to parse from a pileup line
    private static final int expectedTokenCount = 10;
    private static final char fldDelim = '\t';

    // allocate once and don't ever bother creating them again:
    private static final String baseA = "A";
    private static final String baseC = "C";
    private static final String baseG = "G";
    private static final String baseT = "T";
    private static final String emptyStr = ""; // we will use this for "reference" allele in insertions

    /**
     * Return the # of header lines for this file.
     *
     * @param reader the line reader
     * @return 0 in this case, we assume no header lines.
     */
    public Object readHeader(LineReader reader) {
        // we don't require a header line, but it may exist.  We'll deal with that above.
        return null;
    }

    @Override
    public Class<SAMPileupFeature> getFeatureType() {
        return SAMPileupFeature.class;
    }

    public Feature decodeLoc(String line) {
        return decode(line);
    }

    public SAMPileupFeature decode(String line) {
//       0          1             2         3                  4         5            6         7
//*     chrX     466           T         Y                170      170       88      32 ... (piles of read bases  and quals follow)
//*    chrX    141444     *     +CA/+CA       32       468     255     25     +CA     *       5       2       12      6
        String[] tokens = new String[expectedTokenCount];

        // split the line
        int count = ParsingUtils.split(line,tokens,fldDelim);

        // check to see if we've parsed the string into the right number of tokens (expectedTokenCount)
        if (count != expectedTokenCount)
            throw new CodecLineParsingException("the SAM pileup line didn't have the expected number of tokens " +
                                                "(expected = " + expectedTokenCount + ", saw = " + count + " on " +
                                                "line = " + line + ")");

        SAMPileupFeature feature = new SAMPileupFeature();

        feature.setChr(tokens[0]);
        feature.setStart(Integer.parseInt(tokens[1]));

        if(tokens[2].length() != 1)
            throw new CodecLineParsingException("The SAM pileup line had unexpected base " + tokens[2] + " on line = " + line);
        feature.setRef(Character.toUpperCase(tokens[2].charAt(0)));

        String observedString = tokens[3].toUpperCase(); // field 3
        feature.setFWDAlleles(new ArrayList<String>(2));

        feature.setConsensusConfidence(Double.parseDouble(tokens[4]));
        feature.setVariantConfidence(Double.parseDouble(tokens[5]));

        if ( feature.getRef() == '*' ) {
            parseIndels(observedString,feature) ;
            if ( feature.isDeletion() ) feature.setEnd(feature.getStart()+feature.length()-1);
            else feature.setEnd(feature.getStart()); // if it's not a deletion and we are biallelic, this got to be an insertion; otherwise the state is inconsistent!!!!
        } else {
            parseBasesAndQuals(feature,tokens[8],tokens[9]);
            // if the variant is a SNP or a reference base (i.e. no variant at all)
            if ( observedString.length() != 1 ) throw new RuntimeException( "point mutation genotype is expected to be represented by a single letter");
            feature.setRefBases(tokens[2].toUpperCase());
            feature.setEnd(feature.getStart());

            char ch = observedString.charAt(0);

            switch ( ch ) {
                case 'A': feature.getFWDAlleles().add(baseA); feature.getFWDAlleles().add(baseA); break;
                case 'C': feature.getFWDAlleles().add(baseC); feature.getFWDAlleles().add(baseC); break;
                case 'G': feature.getFWDAlleles().add(baseG); feature.getFWDAlleles().add(baseG); break;
                case 'T': feature.getFWDAlleles().add(baseT); feature.getFWDAlleles().add(baseT); break;
                case 'M': feature.getFWDAlleles().add(baseA); feature.getFWDAlleles().add(baseC); break;
                case 'R': feature.getFWDAlleles().add(baseA); feature.getFWDAlleles().add(baseG); break;
                case 'W': feature.getFWDAlleles().add(baseA); feature.getFWDAlleles().add(baseT); break;
                case 'S': feature.getFWDAlleles().add(baseC); feature.getFWDAlleles().add(baseG); break;
                case 'Y': feature.getFWDAlleles().add(baseC); feature.getFWDAlleles().add(baseT); break;
                case 'K': feature.getFWDAlleles().add(baseG); feature.getFWDAlleles().add(baseT); break;
            }
            if ( feature.getFWDAlleles().get(0).charAt(0) == feature.getRef() && feature.getFWDAlleles().get(1).charAt(0) == feature.getRef() ) feature.setVariantType(VariantType.NONE);
            else {
                // 	we know that at least one allele is non-ref;
                // if one is ref and the other is non-ref, or if both are non ref but they are the same (i.e.
                // homozygous non-ref), we still have 2 allelic variants at the site (e.g. one ref and one nonref)
                feature.setVariantType(VariantType.SNP);
                if ( feature.getFWDAlleles().get(0).charAt(0) == feature.getRef() ||
                        feature.getFWDAlleles().get(1).charAt(0) == feature.getRef() ||
                        feature.getFWDAlleles().get(0).equals(feature.getFWDAlleles().get(1))
                        ) feature.setNumNonRef(1);
                else feature.setNumNonRef(2); // if both observations differ from ref and they are not equal to one another, then we get multiallelic site...
            }
        }

        return feature;
    }

    private void parseIndels(String genotype,SAMPileupFeature feature) {
        String [] obs = genotype.split("/"); // get observations, now need to tinker with them a bit

        // if reference allele is among the observed alleles, we will need to take special care of it since we do not have direct access to the reference;
        // if we have an insertion, the "reference" allele is going to be empty; if it it is a deletion, we will deduce the "reference allele" bases
        // from what we have recorded for the deletion allele (e.g. "-CAC")
        boolean hasRefAllele = false;

        for ( int i = 0 ; i < obs.length ; i++ ) {
            if ( obs[i].length() == 1 && obs[i].charAt(0) == '*'  ) {
                hasRefAllele = true;
                feature.getFWDAlleles().add(emptyStr);
                continue;
            }

            String varBases = obs[i].toUpperCase();

            switch ( obs[i].charAt(0) )  {
                case '+':
                    if (!feature.isReference() && !feature.isInsertion()) feature.setVariantType(VariantType.INDEL);
                    else feature.setVariantType(VariantType.INSERTION);
                    feature.setRefBases(emptyStr);
                    break;
                case '-' :
                    if (!feature.isReference() && !feature.isDeletion()) feature.setVariantType(VariantType.INDEL);
                    else feature.setVariantType(VariantType.DELETION);
                    feature.setRefBases(varBases); // remember what was deleted, this will be saved as "reference allele"
                    break;
                default: throw new RuntimeException("Can not interpret observed indel allele record: "+genotype);
            }
            feature.getFWDAlleles().add(varBases);
            feature.setLength(obs[i].length()-1); // inconsistent for non-biallelic indels!!
        }
        if ( hasRefAllele ) {
            // we got at least one ref. allele (out of two recorded)
            if (feature.isReference()) { // both top theories are actually ref allele;
                feature.setNumNonRef(0); // no observations of non-reference allele at all
                feature.setRefBases(emptyStr);
            } else {
                feature.setNumNonRef(1); // hasRefAllele = true, so one allele was definitely ref, hence there is only one left
            }
        } else {
            // we observe two non-ref alleles; they better be the same variant, otherwise the site is not bi-allelic and at the moment we
            // fail to set data in a consistent way.
            if ( feature.getFWDAlleles().get(0).equals(feature.getFWDAlleles().get(1))) feature.setNumNonRef(1);
            else feature.setNumNonRef(2);
        }
        // DONE with indels

    }

    private void parseBasesAndQuals(SAMPileupFeature feature, final String bases, final String quals)
    {
        //System.out.printf("%s%n%s%n", bases, quals);

        // needs to convert the base string with it's . and , to the ref base
        StringBuilder baseBuilder = new StringBuilder();
        StringBuilder qualBuilder = new StringBuilder();
        boolean done = false;
        for ( int i = 0, j = 0; i < bases.length() && ! done; i++ ) {
            //System.out.printf("%d %d%n", i, j);
            char c = (char)bases.charAt(i);

            switch ( c ) {
                case '.':   // matches reference
                case ',':   // matches reference
                    baseBuilder.append(feature.getRef());
                    qualBuilder.append(quals.charAt(j++));
                    break;
                case '$':   // end of read
                    break;
                case '*':   // end of indel?
                    j++;
                    break;
                case '^':   // mapping quality
                    i++;
                    break;
                case '+':   // start of indel
                case '-':   // start of indel
                    final Pattern regex = Pattern.compile("([0-9]+).*");             // matches case 1
                    final String rest = bases.substring(i+1);
                    //System.out.printf("sub is %s%n", rest);
                    Matcher match = regex.matcher(rest);
                    if ( ! match.matches() ) {
                        if ( feature.getRef() != '*' )
                            throw new RuntimeException("Bad pileup format: " + bases + " at position " + i);
                        done = true;
                    }
                    else {
                        String g = match.group(1);
                        //System.out.printf("group is %d, match is %s%n", match.groupCount(), g);
                        int l = Integer.parseInt(g);
                        i += l + g.length();    // length of number + that many bases + +/- at the start (included in the next i++)
                        //System.out.printf("remaining is %d => %s%n", l, bases.substring(i+1));
                    }
                    break;
                default:   // non reference base
                    baseBuilder.append(c);
                    qualBuilder.append(quals.charAt(j++));
            }
        }

        feature.setPileupBases(baseBuilder.toString());
        feature.setPileupQuals(qualBuilder.toString());
    }

}
