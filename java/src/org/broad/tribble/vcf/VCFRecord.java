package org.broad.tribble.vcf;


import org.broad.tribble.Feature;
import org.broad.tribble.util.ParsingUtils;

import java.util.*;

/** the basic VCF record type */
public class VCFRecord implements Feature {

    // the reference base
    private String mReferenceBases;
    // our location
    private String mContig;
    private int mPosition;
    // our id
    private String mID;
    // the alternate bases
    private final List<VCFGenotypeEncoding> mAlts = new ArrayList<VCFGenotypeEncoding>();
    // our qual value
    private double mQual;
    // our filter string
    private String mFilterString;
    // our info fields -- use a TreeMap to ensure they can be pulled out in order (so it passes integration tests)
    private final Map<String, String> mInfoFields = new TreeMap<String, String>();

    // our genotype formatting string
    private String mGenotypeFormatString;

    // the vcf header we're associated with
    private VCFHeader vcfHeader = null;

    // our genotype sample fields
    private final List<VCFGenotypeRecord> mGenotypeRecords = new ArrayList<VCFGenotypeRecord>();

    /**
     * given a reference base, a location, and the format string, create a VCF record.
     *
     * @param referenceBases  the reference bases to use
     * @param contig          our contig
     * @param start           the start location
     * @param genotypeFormatString  the format string
     */
    public VCFRecord(String referenceBases, String contig, int start, String genotypeFormatString) {
        setReferenceBase(referenceBases);
        setLocation(contig, start);
        mGenotypeFormatString = genotypeFormatString;
    }

    /**
     * given the values for each of the columns, create a VCF record.
     *
     * @param columnValues    a mapping of header strings to values
     * @param genotypeFormatString    the format string for the genotype records
     * @param genotypeRecords the genotype records
     */
    public VCFRecord(Map<VCFHeader.HEADER_FIELDS, String> columnValues, String genotypeFormatString, List<VCFGenotypeRecord> genotypeRecords) {
        extractFields(columnValues);
        mGenotypeRecords.addAll(genotypeRecords);
        mGenotypeFormatString = genotypeFormatString;
    }

    /**
     * given the values for each of the columns, create a VCF record.
     *
     * @param columnValues    a mapping of header strings to values
     */
    public VCFRecord(Map<VCFHeader.HEADER_FIELDS, String> columnValues) {
        extractFields(columnValues);
        mGenotypeFormatString = "";
    }

    /**
     * create a VCF record
     *
     * @param referenceBases  the reference bases to use
     * @param contig          the contig this variant is on
     * @param position        our position
     * @param ID              our ID string
     * @param altBases        the list of alternate bases
     * @param qual            the qual field
     * @param filters         the filters used on this variant
     * @param infoFields      the information fields
     * @param genotypeFormatString  the format string
     * @param genotypeObjects the genotype objects
     */
    public VCFRecord(String referenceBases,
                     String contig,
                     long position,
                     String ID,
                     List<VCFGenotypeEncoding> altBases,
                     double qual,
                     String filters,
                     Map<String, String> infoFields,
                     String genotypeFormatString,
                     List<VCFGenotypeRecord> genotypeObjects) {
        setReferenceBase(referenceBases);
        setLocation(contig, position);
        this.mID = ID;
        for (VCFGenotypeEncoding alt : altBases)
            this.addAlternateBase(alt);
        this.setQual(qual);
        this.setFilterString(filters);
        this.mInfoFields.putAll(infoFields);
        this.mGenotypeFormatString = genotypeFormatString;
        this.mGenotypeRecords.addAll(genotypeObjects);
    }

    /**
     * extract the field values from the passed in array
     *
     * @param columnValues a map of the header fields to values
     */
    private void extractFields(Map<VCFHeader.HEADER_FIELDS, String> columnValues) {
        String chrom = null;
        long position = -1;

        for (VCFHeader.HEADER_FIELDS val : columnValues.keySet()) {
            switch (val) {
                case CHROM:
                    chrom = columnValues.get(val);
                    break;
                case POS:
                    position = Integer.valueOf(columnValues.get(val));
                    break;
                case ID:
                    setID(columnValues.get(val));
                    break;
                case REF:
                    if (columnValues.get(val).length() != 1)
                        throw new IllegalArgumentException("Reference base should be a single character");
                    setReferenceBase(columnValues.get(val));
                    break;
                case ALT:
                    String values[] = columnValues.get(val).split(",");
                    for (String alt : values)
                        addAlternateBase(new VCFGenotypeEncoding(alt));
                    break;
                case QUAL:
                    setQual(Double.valueOf(columnValues.get(val)));
                    break;
                case FILTER:
                    setFilterString(columnValues.get(val));
                    break;
                case INFO:
                    String vals[] = columnValues.get(val).split(";");
                    for (String alt : vals) {
                        if ( alt.equals(VCFConstants.EMPTY_INFO_FIELD) )
                            continue;
                        String keyVal[] = alt.split("=");
                        if ( keyVal.length == 1 )
                            addInfoField(keyVal[0], "");
                        else if (keyVal.length == 2)
                            addInfoField(keyVal[0], keyVal[1]);
                        else
                            throw new IllegalArgumentException("info field key-value pair did not parse into key->value pair: " + alt);
                    }
                    break;
            }
        }
        setLocation(chrom, position);
    }

    /**
     * do we have genotyping data
     *
     * @return true if we have genotyping data, false otherwise
     */

    public boolean hasGenotypeData() {
        return (mGenotypeRecords.size() > 0);
    }

    /**
     * @return the ID value for this record
     */
    public String getID() {
        return mID == null ? VCFConstants.EMPTY_ID_FIELD : mID;
    }

    /**
     * get the reference base
     *
     * @return either A, T, C, G, or N
     */
    public String getReference() {
        return mReferenceBases;
    }

    /**
     * get the alternate allele strings
     *
     * @return an array of strings representing the alt alleles, or null if there are none
     */
    public List<String> getAlternateAlleleList() {
        ArrayList<String> alts = new ArrayList<String>();
        for ( VCFGenotypeEncoding alt : mAlts )
            alts.add(alt.getBases());
        return alts;
    }

    public List<VCFGenotypeEncoding> getAlternateAlleles() {
        return mAlts;
    }

    public boolean hasAlternateAllele() {
        for ( VCFGenotypeEncoding alt : mAlts ) {
            if ( alt.getType() != VCFGenotypeEncoding.TYPE.UNCALLED )
                return true;
        }

        return false;
    }

    public boolean isBiallelic() {
        return getAlternateAlleles().size() == 1;
    }

    public boolean isReference() {
        return !hasAlternateAllele();
    }

    public List<String> getAlleleList() {
        ArrayList<String> list = new ArrayList<String>();
        list.add(getReference());
        list.addAll(getAlternateAlleleList());
        return list;
    }

    public double getNonRefAlleleFrequency() {
        if ( mInfoFields.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY) ) {
            return Double.valueOf(mInfoFields.get(VCFConstants.ALLELE_FREQUENCY_KEY));
        } else {
            // this is the poor man's AF
            if ( mInfoFields.containsKey(VCFConstants.ALLELE_COUNT_KEY) && mInfoFields.containsKey(VCFConstants.ALLELE_NUMBER_KEY)) {
                String splt[] = mInfoFields.get(VCFConstants.ALLELE_COUNT_KEY).split(",");
                if ( splt.length > 0 ) {
                    return (Double.valueOf(splt[0]) / Double.valueOf(mInfoFields.get(VCFConstants.ALLELE_NUMBER_KEY)));
                }
            }
        }

        return 0.0;
    }

    public VCFGenotypeEncoding.TYPE getType() {
        VCFGenotypeEncoding.TYPE type = mAlts.get(0).getType();
        for (int i = 1; i < mAlts.size(); i++) {
            if ( mAlts.get(i).getType() != type )
                return VCFGenotypeEncoding.TYPE.MIXED;  // if we have more than one type, return mixed
        }
        return type;
    }

    public boolean isDeletion() {
        return getType() == VCFGenotypeEncoding.TYPE.DELETION;
    }

    public boolean isInsertion() {
        return getType() == VCFGenotypeEncoding.TYPE.INSERTION;
    }

    public boolean isIndel() {
        return isDeletion() || isInsertion();
    }

    public boolean isSNP() {
        return getType() == VCFGenotypeEncoding.TYPE.SINGLE_BASE;
    }

    public boolean isNovel() {
        return ( ! isInDBSNP() ) && ( ! isInHapmap() );
    }

    public boolean isInDBSNP() {
        return ( ( mID != null && ! mID.equals(".") ) || ( mInfoFields.get(VCFConstants.DBSNP_KEY) != null && mInfoFields.get(VCFConstants.DBSNP_KEY).equals("1") ) );
    }

    public boolean isInHapmap() {
        if ( mInfoFields.get(VCFConstants.HAPMAP2_KEY) != null && mInfoFields.get(VCFConstants.HAPMAP2_KEY).equals("1") ) {
            return true;
        } else {
            return ( mInfoFields.get(VCFConstants.HAPMAP3_KEY) != null && mInfoFields.get(VCFConstants.HAPMAP3_KEY).equals("1") );
        }
    }

    public char getAlternativeBaseForSNP() {
        if ( !isSNP() && !isBiallelic() )
            throw new IllegalStateException("This record does not represent a SNP");
        return mAlts.get(0).getBases().charAt(0);
    }

    public char getReferenceForSNP() {
        if ( !isSNP() )
            throw new IllegalStateException("This record does not represent a SNP");
        return getReference().charAt(0);
    }

    /**
     * @return the phred-scaled quality score
     */
    public double getQual() {
        return mQual;
    }

    public int getPosition() {
        return mPosition;
    }

    public boolean isMissingQual() {
        return VCFConstants.MISSING_GENOTYPE_QUALITY_v3.equals(String.valueOf((int)mQual));
    }

    /**
     * @return the -log10PError
     */
    public double getNegLog10PError() {
        return mQual / 10.0;
    }

    /**
     * get the filter criteria
     *
     * @return an array of strings representing the filtering criteria, or UNFILTERED if none are applied
     */
    public String[] getFilteringCodes() {
        if (mFilterString == null) return new String[]{VCFConstants.UNFILTERED};
        return mFilterString.split(VCFConstants.FILTER_CODE_SEPARATOR);
    }

    public boolean isFiltered() {
        String[] codes = getFilteringCodes();
        return !codes[0].equals(VCFConstants.UNFILTERED) && !codes[0].equals(VCFConstants.PASSES_FILTERS_v3);
    }

//    public boolean hasFilteringCodes() {
//        return mFilterString != null;
//    }

    public String getFilterString() {
        return mFilterString;
    }

    /**
     * get the information key-value pairs as a Map<>
     *
     * @return a map, of the info key-value pairs
     */
    public final Map<String, String> getInfoValues() {
        return mInfoFields;
    }

    public List<VCFGenotypeRecord> getVCFGenotypeRecords() {
        return mGenotypeRecords;
    }

    /**
     * @return a List of the sample names
     */
    public String[] getSampleNames() {
        String names[] = new String[mGenotypeRecords.size()];
        for (int i = 0; i < mGenotypeRecords.size(); i++) {
            names[i] = mGenotypeRecords.get(i).getSampleName();
        }
        return names;
    }

   public VCFGenotypeRecord getGenotype(final String sampleName) {
       for ( VCFGenotypeRecord rec : getVCFGenotypeRecords() ) {
           if ( rec.getSampleName().equals(sampleName) ) {
               return rec;
           }
       }

       return null;
   }

    public String getGenotypeFormatString() {
        return mGenotypeFormatString;
    }// the formatting string for our genotype records

    public void setGenotypeFormatString(String newFormatString) {
        mGenotypeFormatString = newFormatString;
    }

    public void setReferenceBase(String reference) {
        mReferenceBases = reference.toUpperCase();
    }

    public void setLocation(String chrom, long position) {
        if ( chrom == null )
            throw new IllegalArgumentException("Chromosomes cannot be missing");
        if ( position < 0 )
            throw new IllegalArgumentException("Position values must be greater than 0");
        this.mContig = chrom;
        this.mPosition = (int)position;
    }

    public void setID(String ID) {
        mID = ID;
    }

    public void setQual(double qual) {
        if ( qual < 0 && !VCFConstants.MISSING_GENOTYPE_QUALITY_v3.equals(String.valueOf((int)qual)) )
            throw new IllegalArgumentException("Qual values cannot be negative unless they are " + VCFConstants.MISSING_GENOTYPE_QUALITY_v3 + " ('unknown')");
        mQual = qual;
    }

    public void setFilterString(String filterString) {
        mFilterString = filterString;
    }

    public void addGenotypeRecord(VCFGenotypeRecord mGenotypeRecord) {
        mGenotypeRecords.add(mGenotypeRecord);
    }

    public void setGenotypeRecords(List<VCFGenotypeRecord> records) {
        mGenotypeRecords.clear();
        for ( VCFGenotypeRecord g : records )
            addGenotypeRecord(g);
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

    public void setAlternateBases(List<VCFGenotypeEncoding> bases) {
        mAlts.clear();
        for ( VCFGenotypeEncoding e : bases )
            addAlternateBase(e);
    }

    /**
     * add an info field to the record
     *
     * @param key   the key, from the spec or a user created key
     * @param value it's value as a string
     */
    public void addInfoField(String key, String value) {
        //System.out.printf("Adding info field %s=%s%n", key, value);
        mInfoFields.put(key, value);
    }

    public void printInfoFields() {
        for ( Map.Entry<String, String> e : mInfoFields.entrySet() ) {
            System.out.printf("  Current info field %s=%s this=%s%n", e.getKey(), e.getValue(), this);
        }
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
     * @param header                the VCF header for this VCF Record
     * @return a string
     */
    public String toStringEncoding(VCFHeader header) {
        StringBuilder builder = new StringBuilder();

        // CHROM \t POS \t ID \t REF \t ALT \t QUAL \t FILTER \t INFO
        builder.append(mContig);
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(mPosition);
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(getID());
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(getReference());
        builder.append(VCFConstants.FIELD_SEPARATOR);
        List<VCFGenotypeEncoding> alts = getAlternateAlleles();
        if ( alts.size() > 0 ) {
            builder.append(alts.get(0));
            for ( int i = 1; i < alts.size(); i++ ) {
                builder.append(",");
                builder.append(alts.get(i));
            }
        } else {
            builder.append(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
        }
        builder.append(VCFConstants.FIELD_SEPARATOR);
        if ( isMissingQual() )
            builder.append(VCFConstants.MISSING_GENOTYPE_QUALITY_v3);
        else
            builder.append(String.format(VCFConstants.DOUBLE_PRECISION_FORMAT_STRING, mQual));
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(ParsingUtils.join(VCFConstants.FILTER_CODE_SEPARATOR, getFilteringCodes()));
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(createInfoString());

        if ( mGenotypeFormatString != null && mGenotypeFormatString.length() > 0 ) {
//            try {
            addGenotypeData(builder, header);
//            } catch (Exception e) {
//                if ( validationStringency == VCFGenotypeWriter.VALIDATION_STRINGENCY.STRICT )
//                    throw new RuntimeException(e);
//            }
        }

        return builder.toString();
    }

    /**
     * create the info string
     *
     * @return a string representing the infomation fields
     */
    protected String createInfoString() {
        StringBuffer info = new StringBuffer();
        boolean isFirst = true;
        for (Map.Entry<String, String> entry : mInfoFields.entrySet()) {
            if ( isFirst )
                isFirst = false;
            else
                info.append(VCFConstants.INFO_FIELD_SEPARATOR);
            info.append(entry.getKey());
            if ( entry.getValue() != null && !entry.getValue().equals("") ) {
                info.append("=");
                info.append(entry.getValue());
            }
        }
        return info.length() == 0 ? VCFConstants.EMPTY_INFO_FIELD : info.toString();
    }

    /**
     * add the genotype data
     *
     * @param builder the string builder
     * @param header  the header object
     */
    private void addGenotypeData(StringBuilder builder, VCFHeader header) {
        Map<String, VCFGenotypeRecord> gMap = genotypeListToMap(getVCFGenotypeRecords());

        StringBuffer tempStr = new StringBuffer();
        if ( header.getGenotypeSamples().size() < getVCFGenotypeRecords().size() ) {
            for ( String sample : gMap.keySet() ) {
                if ( !header.getGenotypeSamples().contains(sample) )
                    System.err.println("Sample " + sample + " is a duplicate or is otherwise not present in the header");
                else
                    header.getGenotypeSamples().remove(sample);
            }
            throw new IllegalStateException("We have more genotype samples than the header specified; please check that samples aren't duplicated");
        }
        tempStr.append(VCFConstants.FIELD_SEPARATOR + mGenotypeFormatString);

        String[] genotypeFormatStrings = mGenotypeFormatString.split(":");

        for ( String genotype : header.getGenotypeSamples() ) {
            tempStr.append(VCFConstants.FIELD_SEPARATOR);
            if ( gMap.containsKey(genotype) ) {
                VCFGenotypeRecord rec = gMap.get(genotype);
                tempStr.append(rec.toStringEncoding(mAlts, genotypeFormatStrings));
                gMap.remove(genotype);
            } else {
                tempStr.append(VCFGenotypeRecord.stringEncodingForEmptyGenotype(genotypeFormatStrings));
            }
        }
        if ( gMap.size() != 0 ) {
            for ( String sample : gMap.keySet() )
                System.err.println("Sample " + sample + " is being genotyped but isn't in the header.");
            throw new IllegalStateException("We failed to use all the genotype samples; there must be an inconsistancy between the header and records");
        }

        builder.append(tempStr);
    }

    /**
     * compare two VCF records
     *
     * @param other the other VCF record
     * @return true if they're equal
     */
    public boolean equals(VCFRecord other) {
        if (!this.mAlts.equals(other.mAlts)) return false;
        if (!this.mReferenceBases.equals(other.mReferenceBases)) return false;
        if (!this.mContig.equals(other.mContig)) return false;
        if (mPosition != other.mPosition) return false;
        if (!this.mID.equals(other.mID)) return false;
        if (this.mQual != other.mQual) return false;
        if ( this.mFilterString == null ) {
            if ( other.mFilterString != null ) return false;
        } else if ( !this.mFilterString.equals(other.mFilterString) ) return false;
        if (!this.mInfoFields.equals(other.mInfoFields)) return false;
        if (!this.mGenotypeRecords.equals(other.mGenotypeRecords)) return false;
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
        for (int i = 0; i < list.size(); i++) {
            VCFGenotypeRecord rec = list.get(i);
            map.put(rec.getSampleName(), rec);
        }
        return map;
    }

    /** Return the features reference sequence name, e.g chromosome or contig */
    public String getChr() {
        return this.mContig;
    }

    /** Return the start position in 1-based coordinates (first base is 1) */
    public int getStart() {
        return this.mPosition;
    }

    /**
     * Return the end position following 1-based fully closed conventions.  The length of a feature is
     * end - start + 1;
     */
    public int getEnd() {
        return this.mPosition;
    }

    /**
     * set the VCF header we're associated with
     * @param header the header
     */
    void setHeader(VCFHeader header) {
        vcfHeader = header;
    }

    /**
     * get the associated header
     * @return the VCF Header
     */
    public VCFHeader getHeader() {
        return vcfHeader;
    }
}
