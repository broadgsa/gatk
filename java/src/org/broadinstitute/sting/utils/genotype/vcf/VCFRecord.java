package org.broadinstitute.sting.utils.genotype.vcf;


import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

/**
 * the basic VCF record type
 */
public class VCFRecord implements Variation, VariantBackedByGenotype {

    // standard info field keys
    public static final String ANCESTRAL_ALLELE_KEY = "AA";
    public static final String ALLELE_COUNT_KEY = "AC";
    public static final String ALLELE_FREQUENCY_KEY = "AF";
    public static final String ALLELE_NUMBER_KEY = "AN";
    public static final String RMS_BASE_QUALITY_KEY = "BQ";
    public static final String DBSNP_KEY = "DB";
    public static final String DEPTH_KEY = "DP";
    public static final String HAPMAP2_KEY = "H2";
    public static final String RMS_MAPPING_QUALITY_KEY = "MQ";
    public static final String SAMPLE_NUMBER_KEY = "NS";
    public static final String STRAND_BIAS_KEY = "SB";

    // commonly used strings that are in the standard
    public static final String FORMAT_FIELD_SEPERATOR = ":";
    public static final String GENOTYPE_FIELD_SEPERATOR = ":";
    public static final String FIELD_SEPERATOR = "\t";
    public static final String FILTER_CODE_SEPERATOR = ";";
    public static final String INFO_FIELD_SEPERATOR = ";";

    // default values
    public static final String UNFILTERED = ".";
    public static final String PASSES_FILTERS = "0";
    public static final String EMPTY_INFO_FIELD = ".";
    public static final String EMPTY_ID_FIELD = ".";
    public static final String EMPTY_ALLELE_FIELD = ".";
    public static final String DOUBLE_PRECISION_FORMAT_STRING = "%.2f";

    // the reference base
    private char mReferenceBase;
    // our location
    private GenomeLoc mLoc;
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

    private final String mGenotypeFormatString;

    // our genotype sample fields
    private final List<Genotype> mGenotypeRecords = new ArrayList<Genotype>();

    /**
     * given a reference base, a location, and the format string, create a VCF record.
     *
     * @param referenceBase   the reference base to use
     * @param location        our genomic location
     * @param genotypeFormatString  the format string
     */
    public VCFRecord(char referenceBase, GenomeLoc location, String genotypeFormatString) {
        setReferenceBase(referenceBase);
        setLocation(location);
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
     * @param referenceBase   the reference base to use
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
        setReferenceBase(referenceBase);
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
        int position = -1;

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
                    setReferenceBase(columnValues.get(val).charAt(0));
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
     * @return this VCF record's location
     */
    public GenomeLoc getLocation() {
        return mLoc;
    }

    /**
     * @return the ID value for this record
     */
    public String getID() {
        return mID == null ? EMPTY_ID_FIELD : mID;
    }

    /**
     * get the reference base
     *
     * @return either A, T, C, G, or N
     */
    public String getReference() {
        return Character.toString(mReferenceBase);
    }

    /**
     * get the reference base
     *
     * @return either A, T, C, G, or N
     */
    public char getReferenceBase() {
        return mReferenceBase;
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
        if ( mInfoFields.containsKey(ALLELE_FREQUENCY_KEY) ) {
            return Double.valueOf(mInfoFields.get(ALLELE_FREQUENCY_KEY));
        } else {
            // this is the poor man's AF
            if ( mInfoFields.containsKey(ALLELE_COUNT_KEY) && mInfoFields.containsKey(ALLELE_NUMBER_KEY)) {
                String splt[] = mInfoFields.get(ALLELE_COUNT_KEY).split(",");
                if ( splt.length > 0 ) {
                    return (Double.valueOf(splt[0]) / Double.valueOf(mInfoFields.get(ALLELE_NUMBER_KEY)));
                }
            }
        }

        return 0.0;
    }
    
    public VARIANT_TYPE getType() {
        if ( !hasAlternateAllele() )
            return VARIANT_TYPE.REFERENCE;

        VCFGenotypeEncoding.TYPE type = mAlts.get(0).getType();
        for (int i = 1; i < mAlts.size(); i++) {
            if ( mAlts.get(i).getType() != type )
                throw new IllegalStateException("The record contains multiple encoding types");
        }

        switch ( type ) {
            case SINGLE_BASE:
                return VARIANT_TYPE.SNP;
            case UNCALLED:
                // If there are no alt alleles, all of the genotypes are reference or no calls, so we're a reference site
                return VARIANT_TYPE.REFERENCE;
            case DELETION:
                return VARIANT_TYPE.DELETION;
            case INSERTION:
                return VARIANT_TYPE.INSERTION;
        }

        throw new IllegalStateException("The record contains unknown genotype encodings");
    }

    public boolean isDeletion() {
        return getType() == VARIANT_TYPE.DELETION;
    }

    public boolean isInsertion() {
        return getType() == VARIANT_TYPE.INSERTION;
    }

    public boolean isIndel() {
        return isDeletion() || isInsertion();
    }

    public boolean isSNP() {
        return getType() == VARIANT_TYPE.SNP;
    }

    public char getAlternativeBaseForSNP() {
        if ( !isSNP() && !isBiallelic() )
            throw new IllegalStateException("This record does not represent a SNP");
        return mAlts.get(0).getBases().charAt(0);
    }

    public char getReferenceForSNP() {
        if ( !isSNP() )
            throw new IllegalStateException("This record does not represent a SNP");
        return getReferenceBase();
    }

    /**
     * @return the phred-scaled quality score
     */
    public double getQual() {
        return mQual;
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
        if (mFilterString == null) return new String[]{UNFILTERED};
        return mFilterString.split(FILTER_CODE_SEPERATOR);
    }

    public boolean isFiltered() {
        String[] codes = getFilteringCodes();
        return !codes[0].equals(UNFILTERED) && !codes[0].equals(PASSES_FILTERS);
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
    public Map<String, String> getInfoValues() {
        if (mInfoFields.size() < 1) {
            Map<String, String> map = new HashMap<String, String>();
            map.put(".", "");
            return map;
        }
        return mInfoFields;
    }

    public List<VCFGenotypeRecord> getVCFGenotypeRecords() {
        ArrayList<VCFGenotypeRecord> list = new ArrayList<VCFGenotypeRecord>();
        for ( Genotype g : mGenotypeRecords )
            list.add((VCFGenotypeRecord)g);       
        return list;
    }

    public List<Genotype> getGenotypes() {
        return mGenotypeRecords;
    }

    public Genotype getCalledGenotype() {
        if ( mGenotypeRecords == null || mGenotypeRecords.size() != 1 )
            throw new IllegalStateException("There is not one and only one genotype associated with this record");
        VCFGenotypeRecord record = (VCFGenotypeRecord)mGenotypeRecords.get(0);
        if ( record.isEmptyGenotype() )
            return null;
        return record;
    }

    public boolean hasGenotype(DiploidGenotype x) {
        if ( mGenotypeRecords == null )
            return false;
        for ( Genotype g : mGenotypeRecords ) {
            if ( DiploidGenotype.valueOf(g.getBases()).equals(x) )
                return true;
        }
        return false;
    }

    /**
     * @return a List of the sample names
     */
    public String[] getSampleNames() {
        String names[] = new String[mGenotypeRecords.size()];
        for (int i = 0; i < mGenotypeRecords.size(); i++) {
            VCFGenotypeRecord rec = (VCFGenotypeRecord)mGenotypeRecords.get(i);
            names[i] = rec.getSampleName();
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

    public void setReferenceBase(char referenceBase) {
        referenceBase = Character.toUpperCase(referenceBase);
        if (referenceBase != 'A' && referenceBase != 'C' && referenceBase != 'T' && referenceBase != 'G' && referenceBase != 'N')
            throw new IllegalArgumentException("Reference base must be one of A,C,G,T,N, we saw: " + referenceBase);
        mReferenceBase = referenceBase;
    }

    public void setLocation(String chrom, int position) {
        if ( chrom == null )
            throw new IllegalArgumentException("Chromosomes cannot be missing");
        if ( position < 0 )
            throw new IllegalArgumentException("Position values must be greater than 0");
        mLoc = GenomeLocParser.createGenomeLoc(chrom, position);
    }

    public void setLocation(GenomeLoc location) {
        mLoc = location.clone();
    }

    public void setID(String ID) {
        mID = ID;
    }

    public void setQual(double qual) {
        if (qual < 0)
            throw new IllegalArgumentException("Qual values must be greater than 0");
        mQual = qual;
    }

    public void setFilterString(String filterString) {
        mFilterString = filterString;
    }

    public void addGenotypeRecord(VCFGenotypeRecord mGenotypeRecord) {
        mGenotypeRecords.add(mGenotypeRecord);
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
        //System.out.printf("Adding info field %s=%s%n", key, value);
        mInfoFields.put(key, value);

        // remove the empty token if it's present
        if ( mInfoFields.containsKey(".") )
            mInfoFields.remove(".");
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
     * @param header the VCF header for this VCF Record
     * @return a string
     */
    public String toStringEncoding(VCFHeader header) {
        return toStringEncoding(header, VCFGenotypeWriter.VALIDATION_STRINGENCY.STRICT);
    }

    /**
     * the generation of a string representation, which is used by the VCF writer
     *
     * @param header                the VCF header for this VCF Record
     * @param validationStringency  the validation stringency
     * @return a string
     */
    public String toStringEncoding(VCFHeader header, VCFGenotypeWriter.VALIDATION_STRINGENCY validationStringency) {
        StringBuilder builder = new StringBuilder();

        // CHROM \t POS \t ID \t REF \t ALT \t QUAL \t FILTER \t INFO
        builder.append(mLoc.getContig());
        builder.append(FIELD_SEPERATOR);
        builder.append(mLoc.getStart());
        builder.append(FIELD_SEPERATOR);
        builder.append(getID());
        builder.append(FIELD_SEPERATOR);
        builder.append(getReferenceBase());
        builder.append(FIELD_SEPERATOR);
        List<VCFGenotypeEncoding> alts = getAlternateAlleles();
        if ( alts.size() > 0 ) {
            builder.append(alts.get(0));
            for ( int i = 1; i < alts.size(); i++ ) {
                builder.append(",");
                builder.append(alts.get(i));
            }
        } else {
            builder.append(EMPTY_ALLELE_FIELD);
        }
        builder.append(FIELD_SEPERATOR);
        builder.append(String.format(DOUBLE_PRECISION_FORMAT_STRING, getQual()));
        builder.append(FIELD_SEPERATOR);
        builder.append(Utils.join(FILTER_CODE_SEPERATOR, getFilteringCodes()));
        builder.append(FIELD_SEPERATOR);
        builder.append(createInfoString());

        if ( mGenotypeFormatString != null && mGenotypeFormatString.length() > 0 ) {
            try {
                addGenotypeData(builder, header);
            } catch (Exception e) {
                if ( validationStringency == VCFGenotypeWriter.VALIDATION_STRINGENCY.STRICT )
                    throw new RuntimeException(e.getMessage());
            }
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
        for (String str : getInfoValues().keySet()) {
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
        StringBuffer tempStr = new StringBuffer();
        if ( header.getGenotypeSamples().size() < getGenotypes().size() )
            throw new IllegalStateException("We have more genotype samples than the header specified");
        tempStr.append(FIELD_SEPERATOR + mGenotypeFormatString);

        Map<String, VCFGenotypeRecord> gMap = genotypeListToMap(getGenotypes());
        String[] genotypeFormatStrings = mGenotypeFormatString.split(":");

        for ( String genotype : header.getGenotypeSamples() ) {
            tempStr.append(FIELD_SEPERATOR);
            if ( gMap.containsKey(genotype) ) {
                VCFGenotypeRecord rec = gMap.get(genotype);
                tempStr.append(rec.toStringEncoding(mAlts, genotypeFormatStrings));
                gMap.remove(genotype);
            } else {
                tempStr.append(VCFGenotypeRecord.EMPTY_GENOTYPE);
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
        if (this.mReferenceBase != other.mReferenceBase) return false;
        if (!this.mLoc.equals(other.mLoc)) return false;
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
    private static Map<String, VCFGenotypeRecord> genotypeListToMap(List<Genotype> list) {
        Map<String, VCFGenotypeRecord> map = new HashMap<String, VCFGenotypeRecord>();
        for (int i = 0; i < list.size(); i++) {
            VCFGenotypeRecord rec = (VCFGenotypeRecord)list.get(i);
            map.put(rec.getSampleName(), rec);
        }
        return map;
    }

}
