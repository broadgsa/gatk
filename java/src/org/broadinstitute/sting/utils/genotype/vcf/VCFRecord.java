package org.broadinstitute.sting.utils.genotype.vcf;


import org.broadinstitute.sting.utils.Utils;

import java.util.*;

/**
 * the basic VCF record type
 */
public class VCFRecord {
    // commonly used strings that are in the standard
    public static final String FORMAT_FIELD_SEPERATOR = ":";
    public static final String GENOTYPE_FIELD_SEPERATOR = ":";
    public static final String FIELD_SEPERATOR = "\t";
    public static final String FILTER_CODE_SEPERATOR = ";";
    public static final String INFO_FIELD_SEPERATOR = ";";
    public static final String EMPTY_INFO_FIELD = ".";
    public static final String DOUBLE_PRECISION_FORMAT_STRING = "%.2f";
    // the reference base
    private char mReferenceBase;
    // our contig
    private String mChrome;
    // our position
    private int mPosition;
    // our id; set to '.' if not available
    private String mID;
    // the alternate bases
    private final List<VCFGenotypeEncoding> mAlts = new ArrayList<VCFGenotypeEncoding>();
    // our qual value
    private double mQual;
    // our filter string
    private String mFilterString;
    // our info fields
    private final Map<String, String> mInfoFields = new HashMap<String, String>();

    private final String mGenotypeFormatString;

    // our genotype sample fields
    private final List<VCFGenotypeRecord> mGenotypeFields = new ArrayList<VCFGenotypeRecord>();

    /**
     * given a VCF header, and the values for each of the columns, create a VCF record.
     *
     * @param columnValues    a mapping of header strings to values
     * @param formatString    the format string for the genotype records
     * @param genotypeRecords the genotype records
     */
    public VCFRecord(Map<VCFHeader.HEADER_FIELDS, String> columnValues, String formatString, List<VCFGenotypeRecord> genotypeRecords) {
        extractFields(columnValues);
        mGenotypeFields.addAll(genotypeRecords);
        mGenotypeFormatString = formatString;
    }

    /**
     * create a VCF record
     *
     * @param referenceBase   the reference base to use
     * @param contig          the contig this variant is on
     * @param position        our position
     * @param ID              our ID string
     * @param altBases        the list of alternate bases
     * @param qual            the qual field
     * @param filters         the filters used on this variant
     * @param infoFields      the information fields
     * @param genotypeObjects the genotype objects
     */
    public VCFRecord(char referenceBase,
                     String contig,
                     int position,
                     String ID,
                     List<VCFGenotypeEncoding> altBases,
                     double qual,
                     String filters,
                     Map<String, String> infoFields,
                     String genotypeFormatString,
                     List<VCFGenotypeRecord> genotypeObjects) {
        this.setReferenceBase(referenceBase);
        this.mChrome = contig;
        this.setPosition(position);
        this.mID = ID;
        for (VCFGenotypeEncoding alt : altBases)
            this.addAlternateBase(alt);
        this.setQual(qual);
        this.setFilterString(filters);
        this.mInfoFields.putAll(infoFields);
        mGenotypeFormatString = genotypeFormatString;
        this.mGenotypeFields.addAll(genotypeObjects);
    }

    /**
     * given a VCF header, and the values for each of the columns, create a VCF record.
     *
     * @param columnValues a mapping of header strings to values
     */
    public VCFRecord(Map<VCFHeader.HEADER_FIELDS, String> columnValues) {
        extractFields(columnValues);
        mGenotypeFormatString = null;
    }

    /**
     * extract the field values from the passed in array
     *
     * @param columnValues a map of the header fields to values
     */
    private void extractFields(Map<VCFHeader.HEADER_FIELDS, String> columnValues) {
        for (VCFHeader.HEADER_FIELDS val : columnValues.keySet()) {
            switch (val) {
                case CHROM:
                    this.setChomosome(columnValues.get(val));
                    break;
                case POS:
                    this.setPosition(Integer.valueOf(columnValues.get(val)));
                    break;
                case ID:
                    this.setID(columnValues.get(val));
                    break;
                case REF:
                    if (columnValues.get(val).length() != 1)
                        throw new IllegalArgumentException("Reference base should be a single character");
                    this.setReferenceBase(columnValues.get(val).charAt(0));
                    break;
                case ALT:
                    String values[] = columnValues.get(val).split(",");
                    for (String alt : values)
                        addAlternateBase(new VCFGenotypeEncoding(alt));
                    break;
                case QUAL:
                    this.setQual(Double.valueOf(columnValues.get(val)));
                    break;
                case FILTER:
                    this.setFilterString(columnValues.get(val));
                    break;
                case INFO:
                    String vals[] = columnValues.get(val).split(";");
                    for (String alt : vals) {
                        String keyVal[] = alt.split("=");
                        if (keyVal.length == 1 && keyVal[0].equals(".")) {
                            this.addInfoField(keyVal[0], "");
                            break;
                        }
                        if (keyVal.length != 2)
                            throw new IllegalArgumentException("info field key-value pair did not parse into key->value pair: " + alt);
                        this.addInfoField(keyVal[0], keyVal[1]);
                    }
                    break;
            }
        }
    }

    /**
     * do we have genotyping data
     *
     * @return true if we have genotyping data, false otherwise
     */

    public boolean hasGenotypeData() {
        return (mGenotypeFields.size() > 0);
    }

    /**
     * @return the string for the chromosome that this VCF record is associated with
     */
    public String getChromosome() {
        return this.mChrome;
    }


    /**
     * @return this VCF records position on the specified chromosome
     */
    public long getPosition() {
        return this.mPosition;
    }

    /**
     * @return the ID value for this record
     */
    public String getID() {
        return this.mID;
    }

    /**
     * get the reference base
     *
     * @return either A, T, C, G, or N
     */
    public char getReferenceBase() {
        return this.mReferenceBase;
    }

    /**
     * get the alternate allele strings
     *
     * @return an array of strings representing the alt alleles, or null if there are none
     */
    public List<VCFGenotypeEncoding> getAlternateAlleles() {
        return this.mAlts;
    }

    public boolean hasAlternateAllele() {
        return getAlternateAlleles().size() > 0;
    }

    /**
     * @return the phred-scaled quality score
     */
    public double getQual() {
        return this.mQual;
    }

    /**
     * get the filter criteria
     *
     * @return an array of strings representing the filtering criteria, or 0 is none are applied
     */
    public String[] getFilteringCodes() {
        if (mFilterString == null) return new String[]{"0"};
        return this.mFilterString.split(FILTER_CODE_SEPERATOR);
    }

    public boolean isFiltered() {
        String[] codes = getFilteringCodes();
        if ( codes.length > 1 ) return true;
        else if ( codes[0].equals(".") || codes[0].equals("0") ) return false;
        else return true;
    }

    public boolean hasFilteringCodes() {
        // todo --- currently always returns true
        return getFilteringCodes() != null;
    }

    public String getFilterString() {
        return mFilterString;
    }

    /**
     * get the information key-value pairs as a Map<>
     *
     * @return a map, of the info key-value pairs
     */
    public Map<String, String> getInfoValues() {
        if (this.mInfoFields.size() < 1) {
            Map<String, String> map = new HashMap<String, String>();
            map.put(".", "");
            return map;
        }
        return this.mInfoFields;
    }

    /**
     * @return the number of columnsof data we're storing
     */
    public int getColumnCount() {
        if (this.hasGenotypeData()) return mGenotypeFields.size() + VCFHeader.HEADER_FIELDS.values().length;
        return VCFHeader.HEADER_FIELDS.values().length;
    }

    /**
     * return the mapping of the format tags to the specified sample's values
     *
     * @return a VCFGenotypeRecord
     */
    public List<VCFGenotypeRecord> getVCFGenotypeRecords() {
        return this.mGenotypeFields;
    }

    /**
     * @return a List of the sample names
     */
    public String[] getSampleNames() {
        String names[] = new String[mGenotypeFields.size()];
        int index = 0;
        for (VCFGenotypeRecord rec : this.mGenotypeFields)
            names[index++] = rec.getSampleName();
        return names;
    }

    public String getGenotypeFormatString() {
        return mGenotypeFormatString;
    }// the formatting string for our genotype records

    public void setReferenceBase(char referenceBase) {
        // upppercase the character
        referenceBase = (char) ((referenceBase > 96) ? referenceBase - 32 : referenceBase);
        if (referenceBase != 'A' && referenceBase != 'C' && referenceBase != 'T' && referenceBase != 'G' && referenceBase != 'N')
            throw new IllegalArgumentException("Reference base must be one of A,C,G,T,N, we saw: " + referenceBase);
        this.mReferenceBase = referenceBase;
    }

    public void setPosition(int mPosition) {
        if (mPosition < 0)
            throw new IllegalArgumentException("Position values must be greater than 0");
        this.mPosition = mPosition;
    }

    public void setChomosome(String mChrome) {
        this.mChrome = mChrome;
    }

    public void setID(String mID) {
        this.mID = mID;
    }

    public void setQual(double mQual) {
        if (mQual < 0)
            throw new IllegalArgumentException("Qual values must be greater than 0");
        this.mQual = mQual;
    }

    public void setFilterString(String mFilterString) {
        this.mFilterString = mFilterString;
    }

    public void addGenotypeField(VCFGenotypeRecord mGenotypeFields) {
        this.mGenotypeFields.add(mGenotypeFields);
    }

    /**
     * add an alternate base to our alternate base list.  All bases are uppercased
     * before being added to the list.
     *
     * @param base the base to add
     */
    public void addAlternateBase(VCFGenotypeEncoding base) {
        if (!mAlts.contains(base)) mAlts.add(base);
    }

    /**
     * add an info field to the record
     *
     * @param key   the key, from the spec or a user created key
     * @param value it's value as a string
     */
    public void addInfoField(String key, String value) {
        this.mInfoFields.put(key, value);
    }


    /**
     * add an info field to the record
     *
     * @param m A map from info keys to info values
     */
    public void addInfoFields(Map<String,String> m) {
        for ( Map.Entry<String, String> e : m.entrySet() )
            addInfoField(e.getKey(), e.getValue());
    }

    /**
     * the generation of a string representation, which is used by the VCF writer
     *
     * @param header the VCF header for this VCF Record
     * @return a string
     */
    public String toStringRepresentation(VCFHeader header) {
        StringBuilder builder = new StringBuilder();

        // CHROM \t POS \t ID \t REF \t ALT \t QUAL \t FILTER \t INFO
        builder.append(getChromosome());
        builder.append(FIELD_SEPERATOR);
        builder.append(getPosition());
        builder.append(FIELD_SEPERATOR);
        builder.append(getID());
        builder.append(FIELD_SEPERATOR);
        builder.append(getReferenceBase());
        builder.append(FIELD_SEPERATOR);
        String alts = "";
        for (VCFGenotypeEncoding str : this.getAlternateAlleles()) alts += str.toString() + ",";
        builder.append((alts.length() > 0) ? alts.substring(0, alts.length() - 1) + FIELD_SEPERATOR : "." + FIELD_SEPERATOR);
        builder.append(String.format(DOUBLE_PRECISION_FORMAT_STRING, getQual()));
        builder.append(FIELD_SEPERATOR);
        builder.append(Utils.join(FILTER_CODE_SEPERATOR, getFilteringCodes()));
        builder.append(FIELD_SEPERATOR);
        builder.append(createInfoString());

        if (this.hasGenotypeData()) {
            addGenotypeData(builder, header);
        }
        return builder.toString();
    }

    /**
     * create the info string
     *
     * @return a string representing the infomation fields
     */
    protected String createInfoString() {
        String info = "";
        for (String str : this.getInfoValues().keySet()) {
            if (str.equals(EMPTY_INFO_FIELD))
                return EMPTY_INFO_FIELD;
            else
                info += str + "=" + getInfoValues().get(str) + INFO_FIELD_SEPERATOR;
        }
        return (info.contains(INFO_FIELD_SEPERATOR)) ? info.substring(0, info.lastIndexOf(INFO_FIELD_SEPERATOR)) : info;
    }

    /**
     * add the genotype data
     *
     * @param builder the string builder
     * @param header  the header object
     */
    private void addGenotypeData(StringBuilder builder, VCFHeader header) {
        builder.append(FIELD_SEPERATOR + this.getGenotypeFormatString());
        if (header.getGenotypeSamples().size() < getVCFGenotypeRecords().size())
            throw new RuntimeException("We have more genotype samples than the header specified");

        Map<String, VCFGenotypeRecord> gMap = genotypeListToMap(getVCFGenotypeRecords());

        for (String genotype : header.getGenotypeSamples()) {
            builder.append(FIELD_SEPERATOR);
            if (gMap.containsKey(genotype)) {
                VCFGenotypeRecord rec = gMap.get(genotype);
                if (!rec.toGenotypeString(this.mAlts).equals(""))
                    builder.append(rec.toGenotypeString(this.mAlts));
                for (String s : rec.getFields().keySet()) {
                    if (rec.getFields().get(s).equals("")) continue;
                    builder.append(GENOTYPE_FIELD_SEPERATOR);
                    builder.append(rec.getFields().get(s));
                }
                gMap.remove(genotype);
            } else {
                builder.append(VCFGenotypeRecord.EMPTY_GENOTYPE);
            }
        }
        if (gMap.size() != 0) {
            for (String sample : gMap.keySet())
                System.err.println("Sample " + sample + " is being genotyped but isn't in the header.");
            throw new RuntimeException("We failed to use all the genotype samples; there must be an inconsistancy between the header and records");
        }
    }

    /**
     * compare two VCF records
     *
     * @param other the other VCF record
     * @return true if they're equal
     */
    public boolean equals(VCFRecord other) {
        if (!this.mAlts.equals(other.mAlts)) return false;
        if (this.mReferenceBase != other.mReferenceBase) return false;
        if (!this.mChrome.equals(other.mChrome)) return false;
        if (this.mPosition != other.mPosition) return false;
        if (!this.mID.equals(other.mID)) return false;
        if (this.mQual != other.mQual) return false;
        if (!this.mFilterString.equals(other.mFilterString)) return false;
        if (!this.mInfoFields.equals(other.mInfoFields)) return false;
        if (!this.mGenotypeFields.equals(other.mGenotypeFields)) return false;
        return true;
    }

    /**
     * create a genotype mapping from a list and their sample names
     *
     * @param list a list of genotype samples
     * @return a mapping of the sample name to VCF genotype record
     */
    private static Map<String, VCFGenotypeRecord> genotypeListToMap(List<VCFGenotypeRecord> list) {
        Map<String, VCFGenotypeRecord> map = new HashMap<String, VCFGenotypeRecord>();
        for (VCFGenotypeRecord rec : list) {
            map.put(rec.getSampleName(), rec);
        }
        return map;
    }

}
