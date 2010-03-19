package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFGenotypeRecord
 *         <p/>
 */
public class VCFGenotypeRecord {

    // key names
    public static final String GENOTYPE_KEY = "GT";
    public static final String GENOTYPE_QUALITY_KEY = "GQ";
    public static final String DEPTH_KEY = "DP";
    public static final String HAPLOTYPE_QUALITY_KEY = "HQ";
    public static final String GENOTYPE_FILTER_KEY = "FT";
    public static final String GENOTYPE_POSTERIORS_TRIPLET_KEY = "GL";
    public static final String OLD_DEPTH_KEY = "RD";

    // the values for empty fields
    public static final String EMPTY_GENOTYPE = "./.";
    public static final String EMPTY_ALLELE = ".";
    public static final int MISSING_GENOTYPE_QUALITY = -1;
    public static final int MISSING_DEPTH = -1;
    public static final int MISSING_HAPLOTYPE_QUALITY = -1;
    public static final String UNFILTERED = ".";

    public static final double MAX_QUAL_VALUE = 99.0;

    // what kind of phasing this genotype has
    public enum PHASE {
        UNPHASED, PHASED, PHASED_SWITCH_PROB, UNKNOWN
    }

    // our record
    private VCFRecord mRecord;

    // our phasing
    private PHASE mPhaseType;

    // our bases(s)
    private final List<VCFGenotypeEncoding> mGenotypeAlleles = new ArrayList<VCFGenotypeEncoding>();

    // our mapping of the format mFields to values
    private final Map<String, String> mFields = new HashMap<String, String>();

    // our sample name
    private String mSampleName;

    /**
     * Create a VCF genotype record
     *
     * @param sampleName  sample name
     * @param genotypes   list of genotypes
     * @param phasing     phasing
     */
    public VCFGenotypeRecord(String sampleName, List<VCFGenotypeEncoding> genotypes, PHASE phasing) {
        mSampleName = sampleName;
        if (genotypes != null)
            this.mGenotypeAlleles.addAll(genotypes);
        mPhaseType = phasing;
    }

    public void setVCFRecord(VCFRecord record) {
        mRecord = record;
    }

    public void setSampleName(String name) {
        mSampleName = name;
    }

    /**
     * Adds a field to the genotype record.
     * Throws an exception if the key is GT, as that's computed internally.
     *
     * @param key    the field name (use static variables above for common fields)
     * @param value  the field value
     */
    public void setField(String key, String value) {
        // make sure the GT field isn't being set
        if ( key.equals(GENOTYPE_KEY) )
            throw new IllegalArgumentException("Setting the GT field is not allowed as that's done internally");
        mFields.put(key, value);
    }

    /**
     * determine the phase of the genotype
     *
     * @param phase the string that contains the phase character
     *
     * @return the phase
     */
    static PHASE determinePhase(String phase) {
        // find the phasing information
        if (phase.equals("/"))
            return PHASE.UNPHASED;
        else if (phase.equals("|"))
            return PHASE.PHASED;
        else if (phase.equals("\\"))
            return PHASE.PHASED_SWITCH_PROB;
        else
            throw new IllegalArgumentException("Unknown genotype phasing parameter");
    }


    public PHASE getPhaseType() {
        return mPhaseType;
    }

    public String getSampleName() {
        return mSampleName;
    }

    public List<VCFGenotypeEncoding> getAlleles() {
        return mGenotypeAlleles;
    }

    public Map<String, String> getFields() {
        return mFields;
    }

    /**
     * @return the phred-scaled quality score
     */
    public double getQual() {
        return ( mFields.containsKey(GENOTYPE_QUALITY_KEY) ? Double.valueOf(mFields.get(GENOTYPE_QUALITY_KEY)) : MISSING_GENOTYPE_QUALITY);
    }

    public boolean isMissingQual() {
        return (int)getQual() == MISSING_GENOTYPE_QUALITY;
    }

    public double getNegLog10PError() {
        double qual = getQual();
        return (qual == MISSING_GENOTYPE_QUALITY ? MISSING_GENOTYPE_QUALITY : qual / 10.0);
    }

    public int getReadCount() {
        return ( mFields.containsKey(DEPTH_KEY) ? Integer.valueOf(mFields.get(DEPTH_KEY)) : MISSING_DEPTH);
    }

    public GenomeLoc getLocation() {
        return mRecord != null ? mRecord.getLocation() : null;
    }

    public String getReference() {
        return mRecord != null ? mRecord.getReference() : "N";
    }

    public String getBases() {
        String genotype = "";
        for ( VCFGenotypeEncoding encoding : mGenotypeAlleles )
            genotype += encoding.getBases();
        return genotype;
    }

    public boolean isVariant(char ref) {
        for ( VCFGenotypeEncoding encoding : mGenotypeAlleles ) {
            if ( encoding.getType() == VCFGenotypeEncoding.TYPE.UNCALLED )
                continue;
            if ( encoding.getType() != VCFGenotypeEncoding.TYPE.SINGLE_BASE ||
                    encoding.getBases().charAt(0) != ref )
                return true;
        }
        return false;
    }

    public boolean isPointGenotype() {
        return (mRecord != null ? !mRecord.isIndel() : true);
    }

    public boolean isHom() {
        if ( mGenotypeAlleles.size() == 0 )
            return true;

        String bases = mGenotypeAlleles.get(0).getBases();
        for ( int i = 1; i < mGenotypeAlleles.size(); i++ ) {
            if ( !bases.equals(mGenotypeAlleles.get(1).getBases()) )
                return false;
        }
        return true;
    }

    public boolean isHet() {
        return !isHom();
    }

    public boolean isNoCall() {
        for ( VCFGenotypeEncoding encoding : mGenotypeAlleles ) {
            if ( encoding.getType() != VCFGenotypeEncoding.TYPE.UNCALLED )
                return false;
        }
        return true;
    }

    public boolean isFiltered() {
        return ( mFields.get(GENOTYPE_FILTER_KEY) != null && ! mFields.get(GENOTYPE_FILTER_KEY).equals("0"));
    }

    public int getPloidy() {
        return 2;
    }

    public VCFRecord getRecord() {
        return mRecord;
    }

    private String toGenotypeString(List<VCFGenotypeEncoding> altAlleles) {
        String str = "";
        boolean first = true;
        for (VCFGenotypeEncoding allele : mGenotypeAlleles) {
            if (allele.getType() == VCFGenotypeEncoding.TYPE.UNCALLED)
                str += VCFGenotypeRecord.EMPTY_ALLELE;
            else
                str += String.valueOf((altAlleles.contains(allele)) ? altAlleles.indexOf(allele) + 1 : 0);
            if (first) {
                switch (mPhaseType) {
                    case UNPHASED:
                        str += "/";
                        break;
                    case PHASED:
                        str += "|";
                        break;
                    case PHASED_SWITCH_PROB:
                        str += "\\";
                        break;
                    case UNKNOWN:
                        throw new UnsupportedOperationException("Unknown phase type");
                }
                first = false;
            }
        }
        return str;

    }

    @Override
    public String toString() {
        return String.format("[VCFGenotype %s %s %s %s]", getLocation(), mSampleName, this.mGenotypeAlleles, mFields);
    }

    public boolean isEmptyGenotype() {
        for ( VCFGenotypeEncoding encoding : mGenotypeAlleles ) {
            if ( encoding.getType() != VCFGenotypeEncoding.TYPE.UNCALLED )
                return false;
        }
        return true;
    }

    public boolean equals(Object other) {
        if (other instanceof VCFGenotypeRecord) {
            if (((VCFGenotypeRecord) other).mPhaseType != this.mPhaseType) return false;
            if (!((VCFGenotypeRecord) other).mGenotypeAlleles.equals(this.mGenotypeAlleles)) return false;
            if (!((VCFGenotypeRecord) other).mFields.equals(mFields)) return false;
            if (!((VCFGenotypeRecord) other).mSampleName.equals(this.mSampleName)) return false;
            return true;
        }
        return false;
    }

    /**
     * output a string representation of the VCFGenotypeRecord, given the alternate alleles
     *
     * @param altAlleles the alternate alleles, needed for toGenotypeString()
     * @param genotypeFormatStrings  genotype format strings
     *
     * @return a string
     */
    public String toStringEncoding(List<VCFGenotypeEncoding> altAlleles, String[] genotypeFormatStrings) {
        StringBuilder builder = new StringBuilder();
        builder.append(toGenotypeString(altAlleles));

        for ( String field : genotypeFormatStrings ) {
            if ( field.equals(GENOTYPE_KEY) )
                continue;

            String value = mFields.get(field);
            if ( value == null && field.equals(OLD_DEPTH_KEY) )
                value = mFields.get(DEPTH_KEY);

            builder.append(VCFRecord.GENOTYPE_FIELD_SEPERATOR);
            if ( value == null || value.equals("") )
                builder.append(getMissingFieldValue(field));
            else
                builder.append(value);
        }

        return builder.toString();
    }

    /**
     * output a string representation of an empty genotype
     *
     * @param genotypeFormatStrings  genotype format strings
     *
     * @return a string
     */
    public static String stringEncodingForEmptyGenotype(String[] genotypeFormatStrings) {
        StringBuilder builder = new StringBuilder();
        builder.append(EMPTY_GENOTYPE);

        for ( String field : genotypeFormatStrings ) {
            if ( field.equals(GENOTYPE_KEY) )
                continue;

            builder.append(VCFRecord.GENOTYPE_FIELD_SEPERATOR);
            builder.append(getMissingFieldValue(field));
        }

        return builder.toString();
    }

    public static String getMissingFieldValue(String field) {
        String result;
        if ( field.equals(GENOTYPE_QUALITY_KEY) )
            result = String.valueOf(MISSING_GENOTYPE_QUALITY);
        else if ( field.equals(DEPTH_KEY) || field.equals(OLD_DEPTH_KEY) )
            result = String.valueOf(MISSING_DEPTH);
        else if ( field.equals(GENOTYPE_FILTER_KEY) )
            result = UNFILTERED;
        else if ( field.equals(GENOTYPE_POSTERIORS_TRIPLET_KEY) )
            result = "0,0,0";
        // TODO -- support haplotype quality
        //else if ( field.equals(HAPLOTYPE_QUALITY_KEY) )
        //    result = String.valueOf(MISSING_HAPLOTYPE_QUALITY);
        else
            result = "-1";
        return result;
    }

    public static Set<VCFFormatHeaderLine> getSupportedHeaderStrings() {
        Set<VCFFormatHeaderLine> result = new HashSet<VCFFormatHeaderLine>();
        result.add(new VCFFormatHeaderLine(GENOTYPE_KEY, 1, VCFFormatHeaderLine.INFO_TYPE.String, "Genotype"));
        result.add(new VCFFormatHeaderLine(GENOTYPE_QUALITY_KEY, 1, VCFFormatHeaderLine.INFO_TYPE.Float, "Genotype Quality"));
        result.add(new VCFFormatHeaderLine(DEPTH_KEY, 1, VCFFormatHeaderLine.INFO_TYPE.Integer, "Read Depth (only filtered reads used for calling)"));
        result.add(new VCFFormatHeaderLine(GENOTYPE_POSTERIORS_TRIPLET_KEY, 3, VCFFormatHeaderLine.INFO_TYPE.Float, "Log-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic"));
        //result.add(new VCFFormatHeaderLine(HAPLOTYPE_QUALITY_KEY, 1, VCFFormatHeaderLine.INFO_TYPE.Integer, "Haplotype Quality"));
        return result;
    }

    public void replaceFields(HashMap<String,String> newFields) {
        mFields.clear();
        for ( String s : newFields.keySet() ) {
            mFields.put(s,newFields.get(s));    
        }
    }
}
