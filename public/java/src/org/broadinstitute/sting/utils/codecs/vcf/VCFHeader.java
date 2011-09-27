package org.broadinstitute.sting.utils.codecs.vcf;


import org.broad.tribble.util.ParsingUtils;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFHeader
 *         <p/>
 *         A class representing the VCF header
 */
public class VCFHeader {

    // the mandatory header fields
    public enum HEADER_FIELDS {
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    }

    // the associated meta data
    private final Set<VCFHeaderLine> mMetaData;
    private final Map<String, VCFInfoHeaderLine> mInfoMetaData = new HashMap<String, VCFInfoHeaderLine>();
    private final Map<String, VCFFormatHeaderLine> mFormatMetaData = new HashMap<String, VCFFormatHeaderLine>();
    private final Map<String, VCFHeaderLine> mOtherMetaData = new HashMap<String, VCFHeaderLine>();

    // the list of auxillary tags
    private final Set<String> mGenotypeSampleNames = new LinkedHashSet<String>();

    // the character string that indicates meta data
    public static final String METADATA_INDICATOR = "##";

    // the header string indicator
    public static final String HEADER_INDICATOR = "#";

    // were the input samples sorted originally (or are we sorting them)?
    private boolean samplesWereAlreadySorted = true;


    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData     the meta data associated with this header
     */
    public VCFHeader(Set<VCFHeaderLine> metaData) {
        mMetaData = new TreeSet<VCFHeaderLine>(metaData);
        loadVCFVersion();
        loadMetaDataMaps();
    }

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData            the meta data associated with this header
     * @param genotypeSampleNames the sample names
     */
    public VCFHeader(Set<VCFHeaderLine> metaData, Set<String> genotypeSampleNames) {
        mMetaData = new TreeSet<VCFHeaderLine>();
        if ( metaData != null )
            mMetaData.addAll(metaData);

        mGenotypeSampleNames.addAll(genotypeSampleNames);

        loadVCFVersion();
        loadMetaDataMaps();

        samplesWereAlreadySorted = ParsingUtils.isSorted(genotypeSampleNames);
    }

    /**
     * Adds a header line to the header metadata.
     *
     * @param headerLine Line to add to the existing metadata component.
     */
    public void addMetaDataLine(VCFHeaderLine headerLine) {
        mMetaData.add(headerLine);    
    }

    /**
     * check our metadata for a VCF version tag, and throw an exception if the version is out of date
     * or the version is not present
     */
    public void loadVCFVersion() {
        List<VCFHeaderLine> toRemove = new ArrayList<VCFHeaderLine>();
        for ( VCFHeaderLine line : mMetaData )
            if ( VCFHeaderVersion.isFormatString(line.getKey())) {
                toRemove.add(line);
            }
        // remove old header lines for now,
        mMetaData.removeAll(toRemove);

    }

    /**
     * load the format/info meta data maps (these are used for quick lookup by key name)
     */
    private void loadMetaDataMaps() {
        for ( VCFHeaderLine line : mMetaData ) {
            if ( line instanceof VCFInfoHeaderLine )  {
                VCFInfoHeaderLine infoLine = (VCFInfoHeaderLine)line;
                mInfoMetaData.put(infoLine.getName(), infoLine);
            }
            else if ( line instanceof VCFFormatHeaderLine ) {
                VCFFormatHeaderLine formatLine = (VCFFormatHeaderLine)line;
                mFormatMetaData.put(formatLine.getName(), formatLine);
            }
            else {
                mOtherMetaData.put(line.getKey(), line);
            }
        }
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
        Set<VCFHeaderLine> lines = new LinkedHashSet<VCFHeaderLine>();
        lines.add(new VCFHeaderLine(VCFHeaderVersion.VCF4_0.getFormatString(), VCFHeaderVersion.VCF4_0.getVersionString()));
        lines.addAll(mMetaData);
        return Collections.unmodifiableSet(lines);
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
        return mGenotypeSampleNames.size() > 0;
    }

    /**
     * were the input samples sorted originally?
     *
     * @return true if the input samples were sorted originally, false otherwise
     */
    public boolean samplesWereAlreadySorted() {
        return samplesWereAlreadySorted;
    }

    /** @return the column count */
    public int getColumnCount() {
        return HEADER_FIELDS.values().length + (hasGenotypingData() ? mGenotypeSampleNames.size() + 1 : 0);
    }

    /**
     * @param key    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFInfoHeaderLine getInfoHeaderLine(String key) {
        return mInfoMetaData.get(key);
    }

    /**
     * @param key    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFFormatHeaderLine getFormatHeaderLine(String key) {
        return mFormatMetaData.get(key);
    }

    /**
     * @param key    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFHeaderLine getOtherHeaderLine(String key) {
        return mOtherMetaData.get(key);
    }
}



