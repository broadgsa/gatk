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

    // our header field ordering, as a linked hash set to guarantee ordering
    private Set<HEADER_FIELDS> mHeaderFields = new LinkedHashSet<HEADER_FIELDS>();

    // the associated meta data
    private final Map<String, String> mMetaData = new HashMap<String, String>();

    // the list of auxillary tags
    private final List<String> auxillaryTags = new ArrayList<String>();

    // the character string that indicates meta data
    public static final String METADATA_INDICATOR = "##";

    // the header string indicator
    public static final String HEADER_INDICATOR = "#";

    /** our log, which we use to capture anything from this class */
    private static Logger logger = Logger.getLogger(VCFHeader.class);

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData
     * @param additionalColumns
     */
    public VCFHeader(Set<HEADER_FIELDS> headerFields, Map<String, String> metaData, List<String> additionalColumns) {
        for (HEADER_FIELDS field : headerFields) mHeaderFields.add(field);
        for (String key : metaData.keySet()) mMetaData.put(key, metaData.get(key));
        for (String col : additionalColumns) auxillaryTags.add(col);
    }

    /**
     * get the header fields in order they're presented in the input file
     *
     * @return a set of the header fields, in order
     */
    public Set<HEADER_FIELDS> getHeaderFields() {
        return mHeaderFields;
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
     * get the auxillary tags
     *
     * @return a list of the extra column names, in order
     */
    public List<String> getAuxillaryTags() {
        return auxillaryTags;
    }
}



