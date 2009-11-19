package org.broadinstitute.sting.utils.genotype.vcf;

import org.apache.log4j.Logger;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFHeader
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public class VCFHeader {

    // the manditory header fields
    public enum HEADER_FIELDS {
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    }

    // the associated meta data
    private final Map<String, String> mMetaData = new HashMap<String, String>();

    // the list of auxillary tags
    private final Set<String> mGenotypeSampleNames = new LinkedHashSet<String>();

    // the character string that indicates meta data
    public static final String METADATA_INDICATOR = "##";

    // the header string indicator
    public static final String HEADER_INDICATOR = "#";

    /** our log, which we use to capture anything from this class */
    private static Logger logger = Logger.getLogger(VCFHeader.class);

    /** do we have genotying data? */
    private boolean hasGenotypingData = false;

    /** the current vcf version we support. */
    private static final String VCF_VERSION = "VCRv3.2";

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData     the meta data associated with this header
     */
    public VCFHeader(Map<String, String> metaData) {
        for (String key : metaData.keySet()) mMetaData.put(key, metaData.get(key));
        checkVCFVersion();
    }

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData            the meta data associated with this header
     * @param genotypeSampleNames the genotype format field, and the sample names
     */
    public VCFHeader(Map<String, String> metaData, Set<String> genotypeSampleNames) {
        for (String key : metaData.keySet()) mMetaData.put(key, metaData.get(key));
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
        if (mMetaData.containsKey("format")) {
            if (mMetaData.get("format").equals(VCF_VERSION))
                return;
            throw new RuntimeException("VCFHeader: VCF version of " + mMetaData.get("format") +
                    " doesn't match the supported version of " + VCF_VERSION);
        }
        throw new RuntimeException("VCFHeader: VCF version isn't present");
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
     * @return a map of the meta data
     */
    public Map<String, String> getMetaData() {
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



