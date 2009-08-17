package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


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
    private final List<Pair<String, String>> mMetaData = new ArrayList<Pair<String, String>>();

    // the list of auxillary tags
    private final List<String> auxillaryTags = new ArrayList<String>();

    // the character string that indicates meta data
    public static final String METADATA_INDICATOR = "##";

    // the header string indicator
    public static final String HEADER_INDICATOR = "#";

    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(VCFHeader.class);

    // patterns we use for detecting meta data and header lines
    private static Pattern pMeta = Pattern.compile("^" + METADATA_INDICATOR + "\\s*(\\S+)\\s*=\\s*(\\S+)\\s*$");


    /**
     * create a VCF header, given an array of strings that all start with at least the # character
     *
     * @param headerStrings a list of header strings
     */
    public VCFHeader(List<String> headerStrings) {
        try {
            Thread.sleep(5000);
        } catch (InterruptedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        // iterate over all the passed in strings
        for (String str : headerStrings) {
            Matcher matcher = pMeta.matcher(str);
            if (matcher.matches()) {
                String metaKey = "";
                String metaValue = "";
                if (matcher.groupCount() < 1) continue;
                if (matcher.groupCount() == 2) metaValue = matcher.group(2);
                metaKey = matcher.group(1);
                mMetaData.add(new Pair<String, String>(metaKey, metaValue));
            }
        }

        // iterate over all the passed in strings
        for (String str : headerStrings) {
            if (str.startsWith("#") && !str.startsWith("##")) {
                String[] strings = str.substring(1).split("\\s+");
                for (String s : strings) {
                    if (mHeaderFields.contains(s)) throw new StingException("Header field duplication is not allowed");
                    try {
                        mHeaderFields.add(HEADER_FIELDS.valueOf(s));
                    } catch (IllegalArgumentException e) {
                        this.auxillaryTags.add(s);
                    }
                }
            }
        }
        if (mHeaderFields.size() != HEADER_FIELDS.values().length) {
            throw new StingException("The VCF header is missing " + (HEADER_FIELDS.values().length - mHeaderFields.size()) + " required fields");
        }
    }

    /**
     * get the header fieldsm in order they're presented in the input file
     * @return
     */
    public Set<HEADER_FIELDS> getHeaderFields() {
        return mHeaderFields;
    }

    /**
     * get the meta data, associated with this header
     * @return
     */
    public List<Pair<String, String>> getMetaData() {
        return mMetaData;
    }

    /**
     * get the auxillary tags
     * @return
     */
    public List<String> getAuxillaryTags() {
        return auxillaryTags;
    }
}



