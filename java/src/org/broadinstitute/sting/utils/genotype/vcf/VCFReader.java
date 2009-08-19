package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

import java.io.*;
import java.nio.charset.Charset;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** The VCFReader class, which given a valid vcf file, parses out the header and VCF records */
public class VCFReader implements Iterator<VCFRecord>, Iterable<VCFRecord> {

    // our VCF header
    private VCFHeader mHeader;

    // our buffered input stream
    private BufferedReader mReader;

    // our next record
    private VCFRecord mNextRecord = null;
    
    // a pattern we use for detecting meta data and header lines
    private static Pattern pMeta = Pattern.compile("^" + VCFHeader.METADATA_INDICATOR + "\\s*(\\S+)\\s*=\\s*(\\S+)\\s*$");

    /**
     * Create a VCF reader, given a VCF file
     *
     * @param vcfFile the vcf file to write
     */
    public VCFReader(File vcfFile) {
        Charset utf8 = Charset.forName("UTF-8");
        try {
            mReader = new BufferedReader(
                    new InputStreamReader(
                            new FileInputStream(vcfFile),
                            utf8));
        } catch (FileNotFoundException e) {
            throw new StingException("VCFReader: Unable to find VCF file: " + vcfFile, e);
        }

        String line = null;
        try {
            ArrayList<String> lines = new ArrayList<String>();
            line = mReader.readLine();
            while (line.startsWith("#")) {
                lines.add(line);
                line = mReader.readLine();
            }
            mHeader = this.createHeader(lines);
            mNextRecord = new VCFRecord(mHeader, line);
        } catch (IOException e) {
            throw new StingException("VCFReader: Failed to parse VCF File on line: " + line, e);
        }

    }

    /** @return true if we have another VCF record to return */
    public boolean hasNext() {
        return (mNextRecord != null);
    }

    /**
     * return the next available VCF record.  Make sure to check availability with a call to hasNext!
     *
     * @return a VCFRecord, representing the next record in the file
     */
    public VCFRecord next() {
        VCFRecord rec = mNextRecord;
        try {
            String line = mReader.readLine();
            if (line == null) mNextRecord = null;
            else mNextRecord = new VCFRecord(mHeader, line);
        } catch (IOException e) {
            mNextRecord = null;
        }
        return rec;
    }

    /** Remove is not supported */
    public void remove() {
        throw new UnsupportedOperationException("Unsupported operation");
    }

    /**
     * create a VCF header, given an array of strings that all start with at least the # character.  This function is
     * package protected so that the VCFReader can access this function
     *
     * @param headerStrings a list of header strings
     */
    protected VCFHeader createHeader(List<String> headerStrings) {

        Map<String,String> metaData = new HashMap<String,String>();
        Set<VCFHeader.HEADER_FIELDS> headerFields = new LinkedHashSet<VCFHeader.HEADER_FIELDS>();
        List<String> auxTags = new ArrayList<String>();
        // iterate over all the passed in strings
        for (String str : headerStrings) {
            Matcher matcher = pMeta.matcher(str);
            if (matcher.matches()) {
                String metaKey = "";
                String metaValue = "";
                if (matcher.groupCount() < 1) continue;
                if (matcher.groupCount() == 2) metaValue = matcher.group(2);
                metaKey = matcher.group(1);
                metaData.put(metaKey, metaValue);
            }
        }

        // iterate over all the passed in strings
        for (String str : headerStrings) {
            if (str.startsWith("#") && !str.startsWith("##")) {
                String[] strings = str.substring(1).split("\\s+");
                for (String s : strings) {
                    if (headerFields.contains(s)) throw new StingException("VCFReader: Header field duplication is not allowed");
                    try {
                        headerFields.add(VCFHeader.HEADER_FIELDS.valueOf(s));
                    } catch (IllegalArgumentException e) {
                        if (!s.equals("FORMAT")) 
                        auxTags.add(s);
                    }
                }
            }
        }
        if (headerFields.size() != VCFHeader.HEADER_FIELDS.values().length) {
            throw new StingException("VCFReader: The VCF column header line is missing " + (VCFHeader.HEADER_FIELDS.values().length - headerFields.size())
                    + " of the " + VCFHeader.HEADER_FIELDS.values().length + " required fields");
        }
        return new VCFHeader(headerFields,metaData,auxTags);
    }

    /**
     *
     * @return get the header associated with this reader
     */
    public VCFHeader getHeader() {
        return this.mHeader;
    }

    @Override
    public Iterator<VCFRecord> iterator() {
        return this;
    }

    public void close() {
        try {
            mReader.close();
        } catch (IOException e) {
            // we don't really care
            Utils.warnUser("Unable to close VCF reader file, this is not fatal, but is worth noting");
        }
    }
}
