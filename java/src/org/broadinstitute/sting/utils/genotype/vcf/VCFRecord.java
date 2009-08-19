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
    private final List<VCFGenotypeRecord> mGenotypeFields;

    // the format String, which specifies what each genotype can contain for values
    private final String mFormatString;

    // the associated header
    private final VCFHeader mHeader;

    /**
     * given a VCF header, and the values for each of the columns, create a VCF record.
     *
     * @param header          the VCF header
     * @param columnValues          a mapping of header strings to values
     * @param formatString    the format string for the genotype records
     * @param genotypeRecords the genotype records
     */
    public VCFRecord(VCFHeader header, Map<VCFHeader.HEADER_FIELDS, String> columnValues, String formatString, List<VCFGenotypeRecord> genotypeRecords) {
        mHeader = header;
        mValues.putAll(columnValues);
        mFormatString = formatString;
        mGenotypeFields = new ArrayList<VCFGenotypeRecord>();
        mGenotypeFields.addAll(genotypeRecords);
    }

    /**
     * given a VCF header, and the values for each of the columns, create a VCF record.
     *
     * @param header          the VCF header
     * @param columnValues          a mapping of header strings to values
     */
    public VCFRecord(VCFHeader header, Map<VCFHeader.HEADER_FIELDS, String> columnValues) {
        mHeader = header;
        mValues.putAll(columnValues);
        mGenotypeFields = null;
        mFormatString = null;
    }

    /**
     * do we have genotyping data
     * @return true if we have genotyping data, false otherwise
     */
    public boolean hasGenotypeData() {
        if (mGenotypeFields==null) {
            return false;
        }
        return true;
    }

    /**
     * get the format string
     * @return the format sting, null if it doesn't exist
     */
    public String getFormatString() {
        return mFormatString;
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
        return mValues.get(VCFHeader.HEADER_FIELDS.CHROM);
    }

    /** @return this VCF records position on the specified chromosome */
    public long getPosition() {
        return Long.valueOf(this.mValues.get(VCFHeader.HEADER_FIELDS.POS));
    }

    /** @return the ID value for this record */
    public String getID() {
        return mValues.get(VCFHeader.HEADER_FIELDS.ID);
    }

    /**
     * get the reference base
     *
     * @return either A, T, C, G, or N
     */
    public char getReferenceBase() {
        // TODO: this field isn't validated correctly
        return mValues.get(VCFHeader.HEADER_FIELDS.REF).charAt(0);
    }

    /**
     * get the alternate allele strings
     *
     * @return an array of strings representing the alt alleles, or null if there are none
     */
    public String[] getAlternateAlleles() {
        if (mValues.get(VCFHeader.HEADER_FIELDS.ALT).trim().equals(".")) {
            return null;
        }
        return mValues.get(VCFHeader.HEADER_FIELDS.ALT).split(",");
    }

    public boolean hasAlternateAllele() {
        return getAlternateAlleles() != null;
    }

    /** @return the phred-scaled quality score */
    public int getQual() {
        return Integer.valueOf(mValues.get(VCFHeader.HEADER_FIELDS.QUAL));
    }

    /**
     * get the filter criteria
     *
     * @return an array of strings representing the filtering criteria, or null if none were applied
     */
    public String[] getFilteringCodes() {
        if (mValues.get(VCFHeader.HEADER_FIELDS.FILTER).trim().equals("0")) {
            return null;
        }
        return mValues.get(VCFHeader.HEADER_FIELDS.ALT).split(";");
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
        return mGenotypeFields.size() + mValues.size();
    }

    /**
     * return the mapping of the format tags to the specified sample's values
     * @return a VCFGenotypeRecord
     */
    public List<VCFGenotypeRecord> getVCFGenotypeRecords() {
        return this.mGenotypeFields;
    }

    /** @return a List of the sample names */
    public String[] getSampleNames() {
        String ret[] = new String[mHeader.getGenotypeSamples().size()];
        mHeader.getGenotypeSamples().toArray(ret);
        return ret;
    }

}
