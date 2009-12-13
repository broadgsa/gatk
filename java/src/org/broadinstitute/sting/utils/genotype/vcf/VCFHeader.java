package org.broadinstitute.sting.utils.genotype.vcf;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFHeader
 *         <p/>
 *         A class representing the VCF header
 */
public class VCFHeader {

    public static final String FILE_FORMAT_KEY = "fileformat";
    public static final String OLD_FILE_FORMAT_KEY = "format";   // from version 3.2


    /** the current vcf version we support. */
    public static final String VCF_VERSION_HEADER = "VCFv";
    public static final String OLD_VCF_VERSION_HEADER = "VCRv";  // from version 3.2
    public static final double VCF_VERSION_NUMBER = 3.3;
    public static final String VCF_VERSION = VCF_VERSION_HEADER + VCF_VERSION_NUMBER;


    // the manditory header fields
    public enum HEADER_FIELDS {
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    }

    // the associated meta data
    private final Set<VCFHeaderLine> mMetaData;

    // the list of auxillary tags
    private final Set<String> mGenotypeSampleNames = new LinkedHashSet<String>();

    // the character string that indicates meta data
    public static final String METADATA_INDICATOR = "##";

    // the header string indicator
    public static final String HEADER_INDICATOR = "#";


    /** do we have genotying data? */
    private boolean hasGenotypingData = false;

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData     the meta data associated with this header
     */
    public VCFHeader(Set<VCFHeaderLine> metaData) {
        mMetaData = new TreeSet<VCFHeaderLine>(metaData);
        checkVCFVersion();
    }

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData            the meta data associated with this header
     * @param genotypeSampleNames the genotype format field, and the sample names
     */
    public VCFHeader(Set<VCFHeaderLine> metaData, Set<String> genotypeSampleNames) {
        mMetaData = new TreeSet<VCFHeaderLine>(metaData);
        for (String col : genotypeSampleNames) {
            if (!col.equals("FORMAT"))
                mGenotypeSampleNames.add(col);
        }
        if (genotypeSampleNames.size() > 0) hasGenotypingData = true;
        checkVCFVersion();
    }

    /**
     * check our metadata for a VCF version tag, and throw an exception if the version is out of date
     * or the version is not present
     */
    public void checkVCFVersion() {
        String version = null;
        for ( VCFHeaderLine line : mMetaData ) {
            if ( line.getKey().equals(FILE_FORMAT_KEY) || line.getKey().equals(OLD_FILE_FORMAT_KEY) ) {
                version = line.getValue();
                break;
            }
        }

        if ( version == null )
            mMetaData.add(new VCFHeaderLine(FILE_FORMAT_KEY, VCF_VERSION));
        else if ( !isSupportedVersion(version) )
            throw new RuntimeException("VCF version " + version +
                    " is not yet supported; only version " + VCF_VERSION + " and earlier can be used");
    }

    private boolean isSupportedVersion(String version) {
        if ( !version.startsWith(VCF_VERSION_HEADER) && !version.startsWith(OLD_VCF_VERSION_HEADER) )
            return false;
        try {
            double dVersion = Double.valueOf(version.substring(VCF_VERSION_HEADER.length()));
            return dVersion <= VCF_VERSION_NUMBER;
        } catch (Exception e) { }

        return false;
    }

    /**
     * get the header fields in order they're presented in the input file (which is now required to be
     * the order presented in the spec).
     *
     * @return a set of the header fields, in order
     */
    public Set<HEADER_FIELDS> getHeaderFields() {
        Set<HEADER_FIELDS> fields = new LinkedHashSet<HEADER_FIELDS>();
        for (HEADER_FIELDS field : HEADER_FIELDS.values())
            fields.add(field);
        return fields;
    }

    /**
     * get the meta data, associated with this header
     *
     * @return a set of the meta data
     */
    public Set<VCFHeaderLine> getMetaData() {
        return mMetaData;
    }

    /**
     * get the genotyping sample names
     *
     * @return a list of the genotype column names, which may be empty if hasGenotypingData() returns false
     */
    public Set<String> getGenotypeSamples() {
        return mGenotypeSampleNames;
    }

    /**
     * do we have genotyping data?
     *
     * @return true if we have genotyping columns, false otherwise
     */
    public boolean hasGenotypingData() {
        return hasGenotypingData;
    }

    /** @return the column count, */
    public int getColumnCount() {
        return HEADER_FIELDS.values().length + ((hasGenotypingData) ? mGenotypeSampleNames.size() + 1 : 0);
    }
}



