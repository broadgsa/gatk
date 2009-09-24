package org.broadinstitute.sting.utils.genotype.vcf;


import org.broadinstitute.sting.utils.Utils;

import java.util.*;

/** the basic VCF record type */
public class VCFRecord {
    public static final String FIELD_SEPERATOR = "\t";
    // the reference base
    private char mReferenceBase;
    // our contig
    private String mChrome;
    // our position
    private int mPosition;
    // our id; set to '.' if not available
    private String mID;
    // the alternate bases
    private final List<String> mAlts = new ArrayList<String>();
    // our qual value
    private int mQual;
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
                     List<String> altBases,
                     int qual,
                     String filters,
                     Map<String, String> infoFields,
                     String genotypeFormatString,
                     List<VCFGenotypeRecord> genotypeObjects) {
        this.setReferenceBase(referenceBase);
        this.mChrome = contig;
        this.setPosition(position);
        this.mID = ID;
        for (String alt : altBases)
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
                        addAlternateBase(alt);
                    break;
                case QUAL:
                    this.setQual(Integer.valueOf(columnValues.get(val)));
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
        if (mGenotypeFields.size() < 1) {
            return false;
        }
        return true;
    }

    /** @return the string for the chromosome that this VCF record is associated with */
    public String getChromosome() {
        return this.mChrome;
    }


    /** @return this VCF records position on the specified chromosome */
    public long getPosition() {
        return this.mPosition;
    }

    /** @return the ID value for this record */
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
    public List<String> getAlternateAlleles() {
        return this.mAlts;
    }

    public boolean hasAlternateAllele() {
        return getAlternateAlleles().size() > 0;
    }

    /** @return the phred-scaled quality score */
    public int getQual() {
        return this.mQual;
    }

    /**
     * get the filter criteria
     *
     * @return an array of strings representing the filtering criteria, or null if none were applied
     */
    public String[] getFilteringCodes() {
        if (mFilterString == null) return new String[]{"0"};
        return this.mFilterString.split(";");
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
        if (this.mInfoFields.size() < 1) {
            Map<String,String> map = new HashMap<String, String>();
            map.put(".","");
            return map;
        }
        return this.mInfoFields;
    }

    /** @return the number of columnsof data we're storing */
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

    /** @return a List of the sample names */
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

    public void setQual(int mQual) {
        if (mQual < 0)
            throw new IllegalArgumentException("Qual values must be greater than 0");
        this.mQual = mQual;
    }

    public void setFilterString(String mFilterString) {
        this.mFilterString = mFilterString;
    }

    public void addGenotypeFields(VCFGenotypeRecord mGenotypeFields) {
        this.mGenotypeFields.add(mGenotypeFields);
    }

    public void addAlternateBase(String base) {
        if (base.length() == 1) {
            char nuc = (char) ((base.charAt(0) > 96) ? base.charAt(0) - 32 : base.charAt(0));
            if (nuc != 'A' && nuc != 'C' && nuc != 'T' && nuc != 'G' && nuc != '.')
                throw new IllegalArgumentException("Alternate base must be either A,C,T,G,. or if an indel it must contain length information: " + base);
        } else {
            // we must be an indel, check that the first character is I or D
            char nuc = (char) ((base.charAt(0) > 96) ? base.charAt(0) - 32 : base.charAt(0));
            if (nuc != 'I' && nuc != 'D')
                throw new IllegalArgumentException("Alternate bases of length greater then one must be an indel: " + base);
        }
        this.mAlts.add(base);
    }

    public void addInfoField(String key, String value) {
        this.mInfoFields.put(key, value);
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();

        // else builder.append(FIELD_SEPERATOR + record.getValue(field));
        // CHROM \t POS \t ID \t REF \t ALT \t QUAL \t FILTER \t INFO
        builder.append(getChromosome() + FIELD_SEPERATOR);
        builder.append(getPosition() + FIELD_SEPERATOR);
        builder.append(getID() + FIELD_SEPERATOR);
        builder.append(getReferenceBase() + FIELD_SEPERATOR);
        String alts = "";
        for (String str : this.getAlternateAlleles()) alts += str + ",";
        builder.append((alts.length() > 0) ? alts.substring(0, alts.length() - 1) + FIELD_SEPERATOR : "." + FIELD_SEPERATOR);
        builder.append(getQual() + FIELD_SEPERATOR);
        builder.append(Utils.join(";", getFilteringCodes()) + FIELD_SEPERATOR);
        String info = "";
        for (String str : this.getInfoValues().keySet()) {
            if (str.equals("."))
                info = ".";
            else
                info += str + "=" + getInfoValues().get(str) + ";";
        }

        if (info.length() > 1) builder.append(info.substring(0, info.length() - 1));
        else builder.append(info);

        if (this.hasGenotypeData()) {
            builder.append(FIELD_SEPERATOR + this.getGenotypeFormatString());
            for (VCFGenotypeRecord rec : this.getVCFGenotypeRecords()) {
                builder.append(FIELD_SEPERATOR);
                if (!rec.toGenotypeString(this.mAlts).equals(""))
                    builder.append(rec.toGenotypeString(this.mAlts));
                for (String s : rec.getFields().keySet()) {
                    if (rec.getFields().get(s).equals("")) continue;
                    builder.append(":");
                    builder.append(rec.getFields().get(s));
                }                
            }
        }
        return builder.toString();
    }

    /**
     * compare two VCF records
     *
     * @param other the other VCF record
     *
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

}
