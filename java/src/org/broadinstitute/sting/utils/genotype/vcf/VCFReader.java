package org.broadinstitute.sting.utils.genotype.vcf;

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

    // our pattern matching for the genotype mFields
    private static final Pattern gtPattern = Pattern.compile("([0-9]+)([\\\\|\\/])([0-9]*)");

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
            throw new RuntimeException("VCFReader: Unable to find VCF file: " + vcfFile, e);
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
            mNextRecord = createRecord(line, mHeader);
        } catch (IOException e) {
            throw new RuntimeException("VCFReader: Failed to parse VCF File on line: " + line, e);
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
            else mNextRecord = createRecord(line, mHeader);
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

        Map<String, String> metaData = new HashMap<String, String>();
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
                    if (headerFields.contains(s))
                        throw new RuntimeException("VCFReader: Header field duplication is not allowed");
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
            throw new RuntimeException("VCFReader: The VCF column header line is missing " + (VCFHeader.HEADER_FIELDS.values().length - headerFields.size())
                    + " of the " + VCFHeader.HEADER_FIELDS.values().length + " required fields");
        }
        return new VCFHeader(headerFields, metaData, auxTags);
    }

    /**
     * create the next VCFRecord, given the input line
     *
     * @param line the line from the file
     *
     * @return the VCFRecord
     */
    public static VCFRecord createRecord(String line, VCFHeader mHeader) {
        // things we need to make a VCF record
        Map<VCFHeader.HEADER_FIELDS, String> values = new HashMap<VCFHeader.HEADER_FIELDS, String>();
        String tokens[] = line.split("\\s+");

        // check to ensure that the column count of tokens is right
        if (tokens.length != mHeader.getColumnCount()) {
            throw new RuntimeException("The input file line doesn't contain enough fields, it should have " + mHeader.getColumnCount() + " fields, it has " + tokens.length);
        }

        int index = 0;
        for (VCFHeader.HEADER_FIELDS field : mHeader.getHeaderFields())
            values.put(field, tokens[index++]);
        // if we have genotyping data, we try and extract the genotype fields
        if (mHeader.hasGenotypingData()) {
            String mFormatString = tokens[index];
            List<VCFGenotypeRecord> genotypeRecords = new ArrayList<VCFGenotypeRecord>();
            index++;
            for (String str : mHeader.getGenotypeSamples()) {
                genotypeRecords.add(getVCFGenotype(str, mFormatString, tokens[index], values.get(VCFHeader.HEADER_FIELDS.ALT).split(","), values.get(VCFHeader.HEADER_FIELDS.REF).charAt(0)));
                index++;
            }
            return new VCFRecord(mHeader, values, mFormatString, genotypeRecords);
        }
        return new VCFRecord(mHeader, values);
    }

    /**
     * generate a VCF genotype record, given it's format string, the genotype string, and allele info
     *
     * @param sampleName     the sample name
     * @param formatString   the format string for this record, which contains the keys for the genotype parameters
     * @param genotypeString contains the phasing information, allele information, and values for genotype parameters
     * @param altAlleles     the alternate allele string array, which we index into based on the field parameters
     * @param referenceBase  the reference base
     */
    public static VCFGenotypeRecord getVCFGenotype(String sampleName, String formatString, String genotypeString, String altAlleles[], char referenceBase) {
        // parameters to create the VCF genotype record
        Map<String,String> tagToValue = new HashMap<String, String>();
        VCFGenotypeRecord.PHASE phase = VCFGenotypeRecord.PHASE.UNKNOWN;
        List<String> bases = new ArrayList<String>();

        String keyStrings[] = formatString.split(":");
        for (String key : keyStrings) {
            String parse;
            int nextDivider;
            if (!genotypeString.contains(":")) {
                nextDivider = genotypeString.length();
                parse = genotypeString;
            } else {
                nextDivider = (genotypeString.indexOf(":") > genotypeString.length()) ? genotypeString.length() : genotypeString.indexOf(":");
                parse = genotypeString.substring(0, nextDivider);
            }
            if (key.equals("GT")) {
                Matcher m = gtPattern.matcher(parse);
                if (!m.matches())
                    throw new RuntimeException("Ubable to match GT genotype flag to it's regular expression");
                phase = VCFGenotypeRecord.determinePhase(m.group(2));
                addAllele(m.group(1),altAlleles,referenceBase,bases);
                if (m.group(3).length() > 0) addAllele(m.group(3),altAlleles,referenceBase,bases);
            }
            tagToValue.put(key,parse);
            if (nextDivider+1 >= genotypeString.length()) nextDivider = genotypeString.length() - 1;
            genotypeString = genotypeString.substring(nextDivider+1,genotypeString.length());
        }
        if (keyStrings.length != tagToValue.size() || genotypeString.length() > 0) throw new RuntimeException("genotype value count doesn't match the key count");
        return new VCFGenotypeRecord(sampleName,tagToValue,bases,phase,referenceBase);
    }


    /**
     * add an alternate allele to the list of alleles we have for a VCF genotype record
     *
     * @param alleleNumber  the allele number, as a string
     * @param altAlleles    the list of alternate alleles
     * @param referenceBase the reference base
     */
    private static void addAllele(String alleleNumber, String[] altAlleles, char referenceBase, List<String> bases) {
        if (Integer.valueOf(alleleNumber) == 0)
            bases.add(String.valueOf(referenceBase));
        else
            bases.add(altAlleles[Integer.valueOf(alleleNumber) - 1]);
    }


    /** @return get the header associated with this reader */
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
