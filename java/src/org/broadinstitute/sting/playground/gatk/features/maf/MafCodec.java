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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.features.maf;

import org.broad.tribble.FeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.LineReader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.IOException;
import java.util.Map;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jan 24, 2011
 * Time: 12:04:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class MafCodec implements FeatureCodec {
     private final static Logger log = Logger.getLogger(MafCodec.class);

     private int expectedTokenCount = -1;

     private int BUILD_COL;
     private int CHR_COL;
     private int START_COL;
     private int END_COL;
     private int REF_ALLELE_COL;
     private int TUMOR_ALLELE1_COL;
     private int TUMOR_ALLELE2_COL;
     private int TUMOR_SAMPLE_COL;
     private int NORMAL_SAMPLE_COL;
    // optional fields (absent from maf lite):
     private int VARTYPE_COL = -1;
     private int STRAND_COL = -1;

     private static String BUILD_COLNAME="NCBI_Build";
     private static String CHR_COLNAME="Chromosome";
     private static String START_COLNAME="Start_position";
     private static String END_COLNAME="End_position";
     private static String REF_ALLELE_COLNAME="Reference_Allele";
     private static String TUMOR_ALLELE1_COLNAME="Tumor_Seq_Allele1";
     private static String TUMOR_ALLELE2_COLNAME="Tumor_Seq_Allele2";
     private static String TUMOR_SAMPLE_COLNAME="Tumor_Sample_Barcode";
     private static String NORMAL_SAMPLE_COLNAME="Matched_Norm_Sample_Barcode";
    // optional fields (absent from maf lite):
     private static String VARTYPE_COLNAME="Variant_Type";
     private static String STRAND_COLNAME="Strand";


     public enum MAF_TYPE {
        UNKNOWN,LITE, ANNOTATED
     }

     private static String INS ="INS";
     private static String DEL ="DEL";
     private static String SNP ="SNP";
     private static String MNP ="MNP";

     private MAF_TYPE mafType=MAF_TYPE.UNKNOWN;

    /**
     * Decode a line to obtain just its FeatureLoc for indexing -- contig, start, and stop.
     *
     * @param line the input line to decode
     * @return  Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public Feature decodeLoc(String line) {
           return decode(line);
    }

    public Feature decode(String line) {

        // ignore commented-out lines
        if (line.startsWith("#")) return null;

        // split the line
        String[] tokens = line.split("\\t");

        if ( expectedTokenCount == -1 ) { // do this only when we receive the first line and do not know the number of columns yet
             // we have not seen a single line yet, let's initialize the number of fields from the first line:
             expectedTokenCount = tokens.length;
             if ( expectedTokenCount == 9 ) {
                mafType = MAF_TYPE.LITE;
                log.info("MAF file appears to be MAF Lite");
             } else {
                 if ( expectedTokenCount == 63 ) {
                     mafType = MAF_TYPE.ANNOTATED;
                     log.info("MAF file appears to be MAF-Annotated");
                 } else {
                     log.info("MAF file has "+expectedTokenCount +" columns, unknown file type");
                 }
             }
             if ( line.contains("Chromosome") && line.contains("Start") && line.contains("Build")) {
                // a naive way to detect the line with column names

                 setColumnsFromHeader(tokens);
                 log.info("MAF file contains header, all required columns found!");
                 return null;
            } else {
                 switch( mafType ) {
                     case UNKNOWN: throw new UserException.MalformedFile("Can not guess type of the MAF file from number of columns and there is no header");
                     case LITE: setMafLiteCols(); break;
                     case ANNOTATED: setMafAnnotatedCols(); break;
                 }
                 log.info("MAF file has no header; assuming standard column order for the MAF type "+mafType);
            }
        }


        if (tokens.length != expectedTokenCount) {
            log.error("MAF line contains wrong number of columns");   
            return null;
        }

        // create a new feature from the line:

        int start = Integer.valueOf(tokens[START_COL]);
        int stop = Integer.valueOf(tokens[END_COL]);

        String eventType="UNKNOWN";

        String ref = tokens[REF_ALLELE_COL];
        String alt1 = tokens[TUMOR_ALLELE1_COL];
        String alt2 = tokens[TUMOR_ALLELE2_COL];

        if ( ref.equals("-") ) {
            // insertion
            eventType = INS;
            stop-- ; // maf lists stop as first base after insertion, convert internally to vcf style
            // perform some format validation:

            if ( alt1.equals("-") && alt2.equals("-") )
                throw new UserException.MalformedFile("Inconsistency in MAF: both alt alleles reported as ref ('-') for an insertion");

            if ( ! alt1.equals("-") && ! alt2.equals("-") && ! alt1.equals(alt2) )
                throw new UserException.MalformedFile("Inconsistency in MAF: two different (non-ref) alt alleles reported for an insertion");

            if ( stop != start )
                throw new UserException.MalformedFile("Inconsistency in MAF: end position for an insertion is not start+1");

        } else {
            if ( alt1.equals("-") || alt2.equals("-") ) {
                // deletion
                eventType = DEL;
                start--; // maf lists start as the first deleted base; convert internally to vcf style
                // perform some format validation:

                if ( ! alt1.equals("-") && ! alt1.equals(ref) )
                    throw new UserException.MalformedFile("Inconsistency in MAF: non-deleted alt allele is not ref for a deletion");

                if ( ! alt2.equals("-") && ! alt2.equals(ref) )
                    throw new UserException.MalformedFile("Inconsistency in MAF: non-deleted alt allele is not ref for a deletion");

                if ( (stop - start) != ref.length() )
                    throw new UserException.MalformedFile("Inconsistency in MAF: deletion length is not end-start+1");

            } else {
                // no '-' alleles --> it's a snp/mnp
                if ( ref.length() == 1 ) {
                    // it's a snp
                    eventType = SNP;
                    if ( stop != start )
                        throw new UserException.MalformedFile("Inconsistency in MAF: start/end positions not equal for a SNP");
                } else {
                    // it's an mnp
                    eventType = MNP;
                    if ( (stop - start + 1) != ref.length() )
                        throw new UserException.MalformedFile("Inconsistency in MAF: MNP length is not end-start+1");
                }

                if ( alt1.length() != ref.length() || alt2.length() != ref.length() )
                    throw new UserException.MalformedFile("Inconsistency in MAF: lengths of ref and alt alleles for a SNP/MNP differ");
                if ( ! alt1.equals(ref) && ! alt2.equals(ref) && ! alt1.equals(alt2) )
                    throw new UserException.MalformedFile("Inconsistency in MAF: two different non-ref alt alleles reported for a SNP/MNP");
            }
        }
        // if we got vartype column, make sure it makes sense:
        if ( VARTYPE_COL != -1 && ! tokens[VARTYPE_COL].equals(eventType) )
            throw new UserException.MalformedFile("Inconsistency in MAF: variant looks like a "+eventType +" but annotated as "+
                tokens[VARTYPE_COL]);

        MafFeature feature = new MafFeature(tokens[CHR_COL],start,stop);
        feature.setVariantType(eventType);
        feature.setRefAllele(ref);
        feature.setObservedTumor(alt1,alt2);
        feature.setTumorSample(tokens[TUMOR_SAMPLE_COL]);
        feature.setNormalSample(tokens[NORMAL_SAMPLE_COL]);

        return feature;
    }

    public Class getFeatureType() {
        return MafFeature.class;
    }


    /**  Read and return the header, or null if there is no header.
     *
     * @return header object
     */
    public Object readHeader(LineReader reader) {
        return null;
    }

    /** Set expected column indices for MafLite
     *
     */
    private void setMafLiteCols() {
        BUILD_COL = 0;
        CHR_COL = 1;
        START_COL = 2;
        END_COL = 3;
        REF_ALLELE_COL = 4;
        TUMOR_ALLELE1_COL = 5;
        TUMOR_ALLELE2_COL = 6;
        TUMOR_SAMPLE_COL = 7;
        NORMAL_SAMPLE_COL = 8;
    }

    private void setMafAnnotatedCols() {
        BUILD_COL = 3;
        CHR_COL = 4;
        START_COL = 5;
        END_COL = 6;
        REF_ALLELE_COL = 10;
        TUMOR_ALLELE1_COL = 11;
        TUMOR_ALLELE2_COL = 12;
        TUMOR_SAMPLE_COL = 15;
        NORMAL_SAMPLE_COL = 16;
        VARTYPE_COL = 9;
        STRAND_COL = 7;
    }

    private void setColumnsFromHeader(String[] tokens) {
        Map<String,Integer> colNames = new HashMap<String,Integer>();
        for ( int i = 0 ; i < tokens.length ; i++ ) colNames.put(tokens[i],i);

        if ( colNames.containsKey(BUILD_COLNAME) ) BUILD_COL = colNames.get(BUILD_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+BUILD_COLNAME+" column");

        if ( colNames.containsKey(CHR_COLNAME) ) CHR_COL = colNames.get(CHR_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+CHR_COLNAME+" column");

        if ( colNames.containsKey(START_COLNAME) ) START_COL = colNames.get(START_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+START_COLNAME+" column");

        if ( colNames.containsKey(END_COLNAME) ) END_COL = colNames.get(END_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+END_COLNAME+" column");

        if ( colNames.containsKey(REF_ALLELE_COLNAME) ) REF_ALLELE_COL = colNames.get(REF_ALLELE_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+REF_ALLELE_COLNAME+" column");

        if ( colNames.containsKey(TUMOR_ALLELE1_COLNAME) ) TUMOR_ALLELE1_COL = colNames.get(TUMOR_ALLELE1_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+TUMOR_ALLELE1_COLNAME+" column");

        if ( colNames.containsKey(TUMOR_ALLELE2_COLNAME) ) TUMOR_ALLELE2_COL = colNames.get(TUMOR_ALLELE2_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+TUMOR_ALLELE2_COLNAME+" column");

        if ( colNames.containsKey(TUMOR_SAMPLE_COLNAME) ) TUMOR_SAMPLE_COL = colNames.get(TUMOR_SAMPLE_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+TUMOR_SAMPLE_COLNAME+" column");

        if ( colNames.containsKey(NORMAL_SAMPLE_COLNAME) ) NORMAL_SAMPLE_COL = colNames.get(NORMAL_SAMPLE_COLNAME);
        else throw new UserException.MalformedFile("Maf file does not have "+NORMAL_SAMPLE_COLNAME+" column");

        // we do not require variant type column but we use it if it's present (for validation):
        if ( colNames.containsKey(VARTYPE_COLNAME) ) VARTYPE_COL = colNames.get(VARTYPE_COLNAME);

        // we do not require strand column but we use it if it's present (for validation):
        if ( colNames.containsKey(STRAND_COLNAME) ) STRAND_COL = colNames.get(STRAND_COLNAME);
    }

}
