package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.Utils;

import java.io.*;
import java.nio.charset.Charset;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

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
        if (vcfFile.getName().endsWith(".gz"))
            openGZipFile(vcfFile);
        else
            openTextVersion(vcfFile);

        String line = null;

        // try and parse the header
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

    /**
     * open a g-zipped version of the VCF format
     *
     * @param vcfGZipFile the file to open
     */
    private void openGZipFile(File vcfGZipFile) {
        try {
            mReader = new BufferedReader(
                    new InputStreamReader(new GZIPInputStream(
                            new FileInputStream(vcfGZipFile))));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("VCFReader: Unable to find VCF file: " + vcfGZipFile, e);
        } catch (IOException e) {
            throw new RuntimeException("VCFReader: A problem occured trying to open the file using the gzipped decompressor, filename: " + vcfGZipFile, e);
        }
    }

    /**
     * open the vcf file as a text file
     *
     * @param vcfFile the vcf file name
     */
    private void openTextVersion(File vcfFile) {
        try {
            mReader = new BufferedReader(
                    new InputStreamReader(
                            new FileInputStream(vcfFile),
                            Charset.forName("UTF-8")));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("VCFReader: Unable to find VCF text file: " + vcfFile, e);
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
     *
     * @return a VCF Header created from the list of stinrgs
     */
    protected VCFHeader createHeader(List<String> headerStrings) {

        Map<String, String> metaData = new HashMap<String, String>();
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
        for (String str : headerStrings) {  // TODO: fix, we shouldn't loop over every line
            if (str.startsWith("#") && !str.startsWith("##")) {
                String[] strings = str.substring(1).split("\\s+");
                // the columns should be in order according to Richard Durbin
                int arrayIndex = 0;
                for (VCFHeader.HEADER_FIELDS field : VCFHeader.HEADER_FIELDS.values()) {
                    try {
                        if (field != VCFHeader.HEADER_FIELDS.valueOf(strings[arrayIndex]))
                            throw new RuntimeException("VCFReader: we were expecting column name " + field + " but we saw " + strings[arrayIndex]);
                    } catch (IllegalArgumentException e) {
                        throw new RuntimeException("VCFReader: Unknown column name \"" + strings[arrayIndex] + "\", it does not match a known column header name.");
                    }
                    arrayIndex++;
                }
                while (arrayIndex < strings.length) {
                    if (!strings[arrayIndex].equals("FORMAT"))
                        auxTags.add(strings[arrayIndex]);
                    arrayIndex++;
                }
            }
        }
        return new VCFHeader(metaData, auxTags);
    }

    /**
     * create the next VCFRecord, given the input line
     *
     * @param line    the line from the file
     * @param mHeader the VCF header
     *
     * @return the VCFRecord
     */
    public static VCFRecord createRecord(String line, VCFHeader mHeader) {
        // things we need to make a VCF record
        Map<VCFHeader.HEADER_FIELDS, String> values = new HashMap<VCFHeader.HEADER_FIELDS, String>();
        String tokens[] = line.split("\\s+");

        // check to ensure that the column count of tokens is right
        if (tokens.length != mHeader.getColumnCount()) {
            throw new RuntimeException("The input file line doesn't contain enough fields, it should have " + mHeader.getColumnCount() + " fields, it has " + tokens.length + ". Line = " + line);
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
     *
     * @return a VCFGenotypeRecord
     */
    public static VCFGenotypeRecord getVCFGenotype(String sampleName, String formatString, String genotypeString, String altAlleles[], char referenceBase) {
        // parameters to create the VCF genotype record
        Map<String, String> tagToValue = new HashMap<String, String>();
        VCFGenotypeRecord.PHASE phase = VCFGenotypeRecord.PHASE.UNKNOWN;
        List<String> bases = new ArrayList<String>();
        int addedCount = 0;
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
                    throw new RuntimeException("VCFReader: Unable to match GT genotype flag to it's expected pattern, the field was: " + parse);
                phase = VCFGenotypeRecord.determinePhase(m.group(2));
                addAllele(m.group(1), altAlleles, referenceBase, bases);
                if (m.group(3).length() > 0) addAllele(m.group(3), altAlleles, referenceBase, bases);
            }
            tagToValue.put(key, parse);
            addedCount++;
            if (nextDivider + 1 >= genotypeString.length()) nextDivider = genotypeString.length() - 1;
            genotypeString = genotypeString.substring(nextDivider + 1, genotypeString.length());
        }
        // catch some common errors, either there are too many field keys or there are two many field values
        if (keyStrings.length != tagToValue.size())
            throw new RuntimeException("VCFReader: genotype value count doesn't match the key count (expected "
                    + keyStrings.length + " but saw " + tagToValue.size() + ")");
        else if (genotypeString.length() > 0)
            throw new RuntimeException("VCFReader: genotype string contained additional unprocessed fields: " + genotypeString
                    + ".  This most likely means that the format string is shorter then the value fields.");
        return new VCFGenotypeRecord(sampleName, tagToValue, bases, phase, referenceBase);
    }


    /**
     * add an alternate allele to the list of alleles we have for a VCF genotype record
     *
     * @param alleleNumber  the allele number, as a string
     * @param altAlleles    the list of alternate alleles
     * @param referenceBase the reference base
     * @param bases         the list of bases for this genotype call
     */
    private static void addAllele(String alleleNumber, String[] altAlleles, char referenceBase, List<String> bases) {
        int alleleValue = Integer.valueOf(alleleNumber);
        // check to make sure the allele value is within bounds
        if (alleleValue < 0 || alleleValue > altAlleles.length)
            throw new IllegalArgumentException("VCFReader: the allele value of " + alleleValue + " is out of bounds given the alternate allele list.");
        if (alleleValue == 0)
            bases.add(String.valueOf(referenceBase));
        else
            bases.add(altAlleles[alleleValue - 1]);
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
