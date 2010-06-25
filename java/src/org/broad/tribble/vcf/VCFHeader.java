package org.broad.tribble.vcf;


import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFHeader
 *         <p/>
 *         A class representing the VCF header
 */
public class VCFHeader {

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

    // our header versionVCF
    private VCFHeaderVersion versionVCF;
    
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
     * check our metadata for a VCF versionVCF tag, and throw an exception if the versionVCF is out of date
     * or the versionVCF is not present
     */
    // TODO: fix this function
    public void checkVCFVersion() {
        VCFHeaderVersion version;
        List<VCFHeaderLine> toRemove = new ArrayList<VCFHeaderLine>();
        for ( VCFHeaderLine line : mMetaData )
            if ( VCFHeaderVersion.isFormatString(line.getKey())) {
                version = VCFHeaderVersion.toHeaderVersion(line.getValue(),line.getKey());
                if (version == null)
                {
                    toRemove.add(line);
                }
                    /**throw new RuntimeException("VCF version " + line.getValue() +
                    " is not supported; only versionVCF " + VCFHeaderVersion.VCF3_2 + " and greater can be used");*/
                else return;
            }
        // remove old header lines for now,
        mMetaData.removeAll(toRemove);
        mMetaData.add(new VCFHeaderLine(VCFHeaderVersion.VCF3_3.getFormatString(), VCFHeaderVersion.VCF3_3.getVersionString()));

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



