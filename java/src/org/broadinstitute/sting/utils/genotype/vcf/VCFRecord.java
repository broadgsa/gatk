package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.StingException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** the basic VCF record type */
public class VCFRecord {
    // required field values
    private final Map<VCFHeader.HEADER_FIELDS, String> mValues = new HashMap<VCFHeader.HEADER_FIELDS, String>();

    // our genotype sample fields
    private final Map<String, String> mGenotypeFields = new HashMap<String, String>();

    // the format String, which specifies what each genotype can contain for values
    private String formatString;

    /**
     * create a VCFRecord, given a VCF header and the the values in this field.  THis is protected, so that the reader is
     * the only accessing object
     *
     * @param header the VCF header
     * @param line   the line to parse into individual fields
     */
    protected VCFRecord(VCFHeader header, String line) {
        String tokens[] = line.split("\\s+");
        List<String> values = new ArrayList<String>();
        for (String str : tokens) values.add(str);
        initialize(header, values);
    }

    /**
     * given a VCF header, and the values for each of the columns, create a VCF record
     *
     * @param header the VCF header
     * @param values the values, as a list, for each of the columns
     */
    public VCFRecord(VCFHeader header, List<String> values) {
        initialize(header, values);
    }

    /**
     * create the VCFRecord
     *
     * @param header the VCF header
     * @param values the list of strings that make up the columns of the record
     */
    private void initialize(VCFHeader header, List<String> values) {
        if (values.size() != header.getColumnCount()) {
            throw new StingException("The input list doesn't contain enough fields, it should have " + header.getColumnCount() + " fields");
        }
        int index = 0;
        for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) {
            mValues.put(field, values.get(index));
            index++;
        }
        if (header.hasGenotypingData()) {
            formatString = values.get(index);
            index++;
            for (String str : header.getGenotypeSamples()) {
                mGenotypeFields.put(str, values.get(index));
                index++;
            }
        }
    }

    /**
     * lookup a value, given it's column name
     *
     * @param key the column name, which is looked up in both the set columns and the auxillary columns
     *
     * @return a String representing the column values, or null if the field doesn't exist in this record
     */
    public String getValue(String key) {
        try {
            return mValues.get(VCFHeader.HEADER_FIELDS.valueOf(key));
        } catch (IllegalArgumentException e) {
            if (this.mGenotypeFields.containsKey(key)) {
                return mGenotypeFields.get(key);
            }
            return null;
        }
    }

    /**
     * get a required field, given the field tag
     *
     * @param field the key for the field
     *
     * @return the field value
     */
    public String getValue(VCFHeader.HEADER_FIELDS field) {
        return mValues.get(field);
    }

    /** @return the string for the chromosome that this VCF record is associated with */
    public String getChromosome() {
        return this.mValues.get(VCFHeader.HEADER_FIELDS.CHROM);
    }

    /** @return this VCF records position on the specified chromosome */
    public long getPosition() {
        return Long.valueOf(this.mValues.get(VCFHeader.HEADER_FIELDS.POS));
    }

    /** @return the ID value for this record */
    public String getID() {
        return this.mValues.get(VCFHeader.HEADER_FIELDS.ID);
    }

    /**
     * get the reference base
     *
     * @return either A, T, C, G, or N
     */
    public char getReferenceBase() {
        // TODO: this field isn't validated correctly
        return this.mValues.get(VCFHeader.HEADER_FIELDS.REF).charAt(0);
    }

    /**
     * get the alternate allele strings
     *
     * @return an array of strings representing the alt alleles, or null if there are none
     */
    public String[] getAlternateAlleles() {
        if (this.mValues.get(VCFHeader.HEADER_FIELDS.ALT).trim().equals(".")) {
            return null;
        }
        return this.mValues.get(VCFHeader.HEADER_FIELDS.ALT).split(",");
    }

    public boolean hasAlternateAllele() {
        return getAlternateAlleles() != null;
    }

    /** @return the phred-scaled quality score */
    public int getQual() {
        return Integer.valueOf(this.mValues.get(VCFHeader.HEADER_FIELDS.QUAL));
    }

    /**
     * get the filter criteria
     *
     * @return an array of strings representing the filtering criteria, or null if none were applied
     */
    public String[] getFilteringCodes() {
        if (this.mValues.get(VCFHeader.HEADER_FIELDS.FILTER).trim().equals("0")) {
            return null;
        }
        return this.mValues.get(VCFHeader.HEADER_FIELDS.ALT).split(";");
    }

    public boolean hasFilteringCodes() {
        return getAlternateAlleles() != null;
    }

    /**
     * get the information key-value pairs as a Map<>
     *
     * @return a map, of the info key-value pairs
     */
    public Map<String, String> getInfoValues() {
        Map<String, String> ret = new HashMap<String, String>();
        String infoSplit[] = mValues.get(VCFHeader.HEADER_FIELDS.INFO).split(";");
        for (String s : infoSplit) {
            String keyValue[] = s.split("=");
            if (keyValue.length != 2)
                throw new StingException("Key value pairs must have both a key and a value; pair: " + s);
            ret.put(keyValue[0], keyValue[1]);
        }
        return ret;
    }

    /** @return the number of columnsof data we're storing */
    public int getColumnCount() {
        return this.mGenotypeFields.size() + this.mValues.size();
    }

    /**
     * return the mapping of the format tags to the specified sample's values
     * @param sampleName the sample name to get the genotyping tags for
     * @return a VCFGenotypeRecord
     */
    public VCFGenotypeRecord getVCFGenotypeRecord(String sampleName) {
        if (!this.mGenotypeFields.containsKey(sampleName)) {
            throw new IllegalArgumentException("Sample Name: " + sampleName + " doesn't exist in this VCF record");    
        }
        return new VCFGenotypeRecord(formatString,mGenotypeFields.get(sampleName),this.getAlternateAlleles(),this.getReferenceBase());

    }

}
