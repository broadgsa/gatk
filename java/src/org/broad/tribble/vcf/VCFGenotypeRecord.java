package org.broad.tribble.vcf;

import org.broadinstitute.sting.utils.Utils;

import java.util.*;


/**
 *
 * @author aaron
 *
 * Class VCFGenotypeRecord
 *
 * the basics of a genotype call in VCF
 */
public class VCFGenotypeRecord {

    public static final double MAX_QUAL_VALUE = 99.0;

    // what kind of phasing this genotype has
    public enum PHASE {
        UNPHASED("/"), PHASED("|"), PHASED_SWITCH_PROB("\\"); // , UNKNOWN

        String genotypeSeparator;
        PHASE(String sep) { this.genotypeSeparator = sep; }
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
        if ( key.equals(VCFConstants.GENOTYPE_KEY) )
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
        for ( PHASE p : PHASE.values() ) {
            if (phase.equals(p.genotypeSeparator))
            return p;
        }

        throw new IllegalArgumentException("Unknown genotype phasing parameter: " + phase);
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
        return ( mFields.containsKey(VCFConstants.GENOTYPE_QUALITY_KEY) ? Double.valueOf(mFields.get(VCFConstants.GENOTYPE_QUALITY_KEY)) : Double.valueOf(VCFConstants.MISSING_GENOTYPE_QUALITY_v3));
    }

    public boolean isMissingQual() {
        return VCFConstants.MISSING_GENOTYPE_QUALITY_v3.equals(String.valueOf((int)getQual()));
    }

    public double getNegLog10PError() {
        return (isMissingQual() ? Double.valueOf(VCFConstants.MISSING_GENOTYPE_QUALITY_v3) : getQual() / 10.0);
    }

    public int getReadCount() {
        return ( mFields.containsKey(VCFConstants.DEPTH_KEY) ? Integer.valueOf(mFields.get(VCFConstants.DEPTH_KEY)) : Integer.valueOf(VCFConstants.MISSING_DEPTH_v3));
    }

    public String getLocation() {
        return mRecord != null ? mRecord.getChr() + ":" + mRecord.getPosition() : null;
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
        return ( mFields.get(VCFConstants.GENOTYPE_FILTER_KEY) != null &&
                !mFields.get(VCFConstants.GENOTYPE_FILTER_KEY).equals(VCFConstants.UNFILTERED) &&
                !mFields.get(VCFConstants.GENOTYPE_FILTER_KEY).equals(VCFConstants.PASSES_FILTERS_v3));
    }

    public int getPloidy() {
        return mGenotypeAlleles.size();
    }

    public VCFRecord getRecord() {
        return mRecord;
    }

    private String toGenotypeString(List<VCFGenotypeEncoding> altAlleles) {
        List<String> alleleStrings = new ArrayList<String>(altAlleles.size());
        for (VCFGenotypeEncoding allele : mGenotypeAlleles) {
            if (allele.getType() == VCFGenotypeEncoding.TYPE.UNCALLED)
                alleleStrings.add(VCFConstants.EMPTY_ALLELE);
            else
                alleleStrings.add(String.valueOf((altAlleles.contains(allele)) ? altAlleles.indexOf(allele) + 1 : 0));
        }

        return Utils.join(mPhaseType.genotypeSeparator, alleleStrings);
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
        return toStringEncoding(altAlleles, genotypeFormatStrings, false);
    }

    public String toStringEncoding(List<VCFGenotypeEncoding> altAlleles, String[] genotypeFormatStrings, boolean doVCF40) {
        StringBuilder builder = new StringBuilder();
        builder.append(toGenotypeString(altAlleles));

        for ( String field : genotypeFormatStrings ) {
            if ( field.equals(VCFConstants.GENOTYPE_KEY) )
                continue;

            String value = mFields.get(field);
            if ( value == null && field.equals(VCFConstants.OLD_DEPTH_KEY) )
                value = mFields.get(VCFConstants.DEPTH_KEY);

            builder.append(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
            if ( value == null || value.equals("") )
                builder.append(getMissingFieldValue(field, doVCF40));
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
        // backward compatibility to VCF 3.3
        return stringEncodingForEmptyGenotype(genotypeFormatStrings, false);
    }
    public static String stringEncodingForEmptyGenotype(String[] genotypeFormatStrings, boolean doVCF40) {
        StringBuilder builder = new StringBuilder();
        builder.append(VCFConstants.EMPTY_GENOTYPE);

        for ( String field : genotypeFormatStrings ) {
            if ( field.equals(VCFConstants.GENOTYPE_KEY) )
                continue;

            // in VCF4.0, if a genotype is empty only the ./. key can be included
            if (!doVCF40) {
                builder.append(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
                builder.append(getMissingFieldValue(field));
            }
        }

        return builder.toString();
    }

    public static String getMissingFieldValue(String field) {
        // backward compatibility to VCF 3.3
        return getMissingFieldValue(field, false);
    }
    public static String getMissingFieldValue(String field, boolean doVCF40) {
        String result;
        if (doVCF40) {
            result = "."; // default missing value
            // TODO - take number of elements in field as input and output corresponding .'s
            if ( field.equals(VCFConstants.GENOTYPE_LIKELIHOODS_KEY) )
                result = ".,.,.";
            else if ( field.equals(VCFConstants.HAPLOTYPE_QUALITY_KEY) )
                result = ".,.";

        }
        else {
            result = "";


            if ( field.equals(VCFConstants.GENOTYPE_QUALITY_KEY) )
                result = String.valueOf(VCFConstants.MISSING_GENOTYPE_QUALITY_v3);
            else if ( field.equals(VCFConstants.DEPTH_KEY) || field.equals(VCFConstants.OLD_DEPTH_KEY) )
                result = String.valueOf(VCFConstants.MISSING_DEPTH_v3);
            else if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY) )
                result = VCFConstants.UNFILTERED;
            else if ( field.equals(VCFConstants.GENOTYPE_LIKELIHOODS_KEY) )
                result = "0,0,0";
            // TODO -- support haplotype quality
            //else if ( field.equals(HAPLOTYPE_QUALITY_KEY) )
            //    result = String.valueOf(MISSING_HAPLOTYPE_QUALITY);
        }
        return result;
    }

    public static Set<VCFFormatHeaderLine> getSupportedHeaderStrings(VCFHeaderVersion version) {
        Set<VCFFormatHeaderLine> result = new HashSet<VCFFormatHeaderLine>();
        result.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
        result.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Float, "Genotype Quality"));
        result.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Read Depth (only filtered reads used for calling)"));
        result.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, 3, VCFHeaderLineType.Float, "Log-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic"));
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
