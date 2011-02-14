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
import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.LineReader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.IOException;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.lang.reflect.Field;

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


     private Column BUILD_COL = new Column("NCBI_Build",true);
     private Column CHR_COL = new Column("Chromosome",true);
     private Column START_COL = new Column("Start_position",true);
     private Column END_COL = new Column("End_position",true);
     private Column REF_ALLELE_COL = new Column("Reference_Allele",true);
     private Column TUMOR_ALLELE1_COL = new Column("Tumor_Seq_Allele1",true);
     private Column TUMOR_ALLELE2_COL = new Column("Tumor_Seq_Allele2",true);
     private Column TUMOR_SAMPLE_COL = new Column("Tumor_Sample_Barcode",true);
     private Column NORMAL_SAMPLE_COL = new Column("Matched_Norm_Sample_Barcode",true);
    // optional fields (absent from maf lite):
     private Column VARTYPE_COL = new Column("Variant_Type",false);
     private Column STRAND_COL = new Column("Strand",false);
     private Column HUGO_GENE_COL = new Column("Hugo_Symbol",false);
     private Column VARCLASS_COL = new Column("Variant_Classification",false);


     public enum MAF_TYPE {
        UNKNOWN,LITE, ANNOTATED
     }

     private static String INS ="INS";
     private static String DEL ="DEL";
     private static String SNP ="SNP";
     private static String MNP ="MNP";

     private MAF_TYPE mafType=MAF_TYPE.UNKNOWN;

     private List<Column> allColumns = null; /// filled dynamically by constructor through introspection. Slow but less typing.

    private boolean tooManyColsWarned = false;
    private boolean tooFewColsWarned = false;

     public MafCodec() {
         allColumns = new ArrayList<Column>(30);
         Field[] fields = this.getClass().getDeclaredFields();
         try {
            for ( Field f : fields ) {
                 if ( f.get(this) instanceof Column ) {
                     allColumns.add((Column)f.get(this));
                 }
            }
         } catch (IllegalAccessException e) {
             throw new StingException("Error in MAFCodec when trying to introspect itself, this is probably a BUG",e);
         }
     }


    /**
     * Decode a line to obtain just its FeatureLoc for indexing -- contig, start, and stop.
     * This method will NOT fill in the additional information available in the maf file
     * @param line the input line to decode
     * @return  Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public Feature decodeLoc(String line) {
           return reallyDecode(line,false);
    }


    /**
     * Fully decode a line, will try extracting as much additional/annotation information from the maf file as it can.
     * @param line the input line to decode
     * @return  Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public Feature decode(String line) {
           return reallyDecode(line,true);
    }

    /** Decodes a maf line. If <code>extra</code> is false, will decode only location and return;
     * if <code>extra</code> is true, then extracts everything it can (samples, annotations, etc)
     * @param line
     * @param extra
     * @return
     */
    public Feature reallyDecode(String line, boolean extra) {

        // ignore commented-out lines
        if (line.startsWith("#")) return null;

        // split the line
        String[] tokens = line.split("\\t",-1);

        if ( expectedTokenCount == -1 ) { // do this only when we receive the first line and do not know the number of columns yet
             // we have not seen a single line yet, let's initialize the number of fields from the first line:
             expectedTokenCount = tokens.length;
             log.info("MAF: line has "+expectedTokenCount+" fields (columns)");
             if ( expectedTokenCount == 9 ) {
                mafType = MAF_TYPE.LITE;
                log.info("MAF file appears to be MAF Lite");
             } else {
                 if ( expectedTokenCount >= 63 ) {
                     mafType = MAF_TYPE.ANNOTATED;
                     log.info("MAF file appears to be MAF-Annotated");
                 } else {
                     log.info("MAF file has "+expectedTokenCount +" columns in first line, unknown file type");
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


        if (tokens.length < expectedTokenCount) {
            if ( ! tooFewColsWarned ) {
                log.error("MAF line contains too few columns ("+tokens.length+"); this error is reported only once.");
                tooFewColsWarned = true;
            }
        }
        if (tokens.length > expectedTokenCount) {
            if ( ! tooManyColsWarned ) {
                log.warn("MAF line contains more columns than expected ("+tokens.length+"); extra columns discarded. This error is shown only once.");
                tooManyColsWarned = true;
            }
        }

        if ( tokens[CHR_COL.getIndex()].equals("Chromosome") ) return null; // if someone uses this codec manually and feeds it the header line multiple times...
        // create a new feature from the line:

        int start = 0;
        try {
            start = Integer.parseInt(START_COL.getValue(tokens));
        } catch (NumberFormatException e) {
            throw new UserException.MalformedFile("Missing or non-numeric start position in line:\n"+line,e);
        }
        int stop = 0 ;
        try {
            stop = Integer.parseInt(END_COL.getValue(tokens));
        } catch (NumberFormatException e) {
            throw new UserException.MalformedFile("Missing or non-numeric stop position in line:\n"+line,e);
        }

        String eventType="UNKNOWN";

        String ref = REF_ALLELE_COL.getValue(tokens);
        String alt1 = TUMOR_ALLELE1_COL.getValue(tokens);
        String alt2 = TUMOR_ALLELE2_COL.getValue(tokens);

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
        if ( VARTYPE_COL.isSet(tokens) && ! tokens[VARTYPE_COL.getIndex()].equals(eventType) )  {
            // special case: we annotate everything as MNP while MAF can have DNP/TNP, these are fine:
            if ( eventType == MNP && (
                    tokens[VARTYPE_COL.getIndex()].equals("DNP") && ref.length() == 2 ||
                    tokens[VARTYPE_COL.getIndex()].equals("TNP") && ref.length() == 3)
                    ) {}                                                              // these are fine
            else {
                throw new UserException.MalformedFile("Inconsistency in MAF: variant looks like a "+eventType +" but annotated as "+
                    tokens[VARTYPE_COL.getIndex()]);
            }
        }
        MafFeature feature = new MafFeature(CHR_COL.getValue(tokens),start,stop);

        if ( ! extra ) return feature; // ignore additional fields unless we were explicitly asked to read those!
        
        feature.setVariantType(eventType);
        feature.setRefAllele(ref);
        feature.setObservedTumor(alt1,alt2);
        feature.setTumorSample(TUMOR_SAMPLE_COL.getValue(tokens));
        feature.setNormalSample(NORMAL_SAMPLE_COL.getValue(tokens));

        if ( HUGO_GENE_COL.isSet(tokens) ) feature.setHugoGeneSymbol(tokens[HUGO_GENE_COL.getIndex()]);
        if ( VARCLASS_COL.isSet(tokens) ) feature.setVariantClassification(tokens[VARCLASS_COL.getIndex()]);

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
        BUILD_COL.setIndex(0);
        CHR_COL.setIndex(1);
        START_COL.setIndex(2);
        END_COL.setIndex(3);
        REF_ALLELE_COL.setIndex(4);
        TUMOR_ALLELE1_COL.setIndex(5);
        TUMOR_ALLELE2_COL.setIndex(6);
        TUMOR_SAMPLE_COL.setIndex(7);
        NORMAL_SAMPLE_COL.setIndex(8);
    }

    private void setMafAnnotatedCols() {
        BUILD_COL.setIndex(3);
        CHR_COL.setIndex(4);
        START_COL.setIndex(5);
        END_COL.setIndex(6);
        REF_ALLELE_COL.setIndex(10);
        TUMOR_ALLELE1_COL.setIndex(11);
        TUMOR_ALLELE2_COL.setIndex(12);
        TUMOR_SAMPLE_COL.setIndex(15);
        NORMAL_SAMPLE_COL.setIndex(16);
        VARTYPE_COL.setIndex(9);
        STRAND_COL.setIndex(7);
        VARCLASS_COL.setIndex(8);
        HUGO_GENE_COL.setIndex(0);
    }

    private void setColumnsFromHeader(String[] tokens) {
        Map<String,Integer> colNames = new HashMap<String,Integer>();
        for ( int i = 0 ; i < tokens.length ; i++ ) colNames.put(tokens[i],i);

        for ( Column c : allColumns ) c.setFromMap(colNames);
    }


}


class Column {
    int index ;
    String name;
    boolean required;

    Column(String name, boolean required) {
        this.name = name;
        this.required = required;
        this.index = -1;
    }

    public String getName() { return name; }
    public void setName(String name) { this.name = name; }
    public int getIndex() { return index; }
    public void setIndex(int index) { this.index = index; }
    public String getValue(String[] fields) {
        if ( index < fields.length ) return fields[index];

        if ( required ) throw new UserException.MalformedFile("In MAF file: required column "+name+" has index "+index+
                    ", but only "+fields.length+ " fields are present in maf line");
        return null;
    }

    /** Sets this column's index from the provided name->index map (i.e. searches for itself in the map).
     * If column not found, <code>throw_exception</code> is true <i>AND</i> this column is required, then an exception will
     * be thrown right away; otherwise returns quietely even if map does not contain this column.
     * @param m
     * @param throw_exception
     */
    public void setFromMap(Map<String,Integer> m, boolean throw_exception) {
        Integer i = m.get(this.name);
        if ( i == null ) {
            if ( this.required && throw_exception ) throw new UserException.MalformedFile("Required column "+this.name+" is missing from the maf file");
            index = -1;
            return; // not found
        }
        this.index = i.intValue(); // found and set.
    }

/**  Sets this column's index from the provided name->index map (i.e. searches for itself in the map).
 * If this column is required but not found in the map, then an exception will
 * be thrown.
 * @param m
 */
    public void setFromMap(Map<String,Integer> m) {
        setFromMap(m,true);
    }

    public boolean isSet() { return index > -1; }

    public boolean isSet(String[] fields) { return index > -1 && index < fields.length; }

}

