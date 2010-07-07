package org.broad.tribble.vcf;



import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** The VCFReaderUtils class, which contains a collection of utilities for working with VCF files */
public class VCFReaderUtils {

    // our pattern matching for the genotype mFields
    private static final Pattern gtPattern = Pattern.compile("([0-9\\.]+)([\\\\|\\/])([0-9\\.]*)");

    /**
     * create a VCF header, given an array of strings that all start with at least the # character.  This function is
     * package protected so that the VCFReaderUtils can access this function
     *
     * @param headerStrings a list of header strings
     *
     * @return a VCF Header created from the list of stinrgs
     */
    public static VCFHeader createHeader(List<String> headerStrings, VCFHeaderVersion version) {
        Set<VCFHeaderLine> metaData = new TreeSet<VCFHeaderLine>();
        Set<String> auxTags = new LinkedHashSet<String>();
        // iterate over all the passed in strings
        for ( String str : headerStrings ) {
            if ( !str.startsWith("##") ) {
                String[] strings = str.substring(1).split("\\t");
                // the columns should be in order according to Richard Durbin
                int arrayIndex = 0;
                for (VCFHeader.HEADER_FIELDS field : VCFHeader.HEADER_FIELDS.values()) {
                    try {
                        if (field != VCFHeader.HEADER_FIELDS.valueOf(strings[arrayIndex]))
                            throw new RuntimeException("VCFReaderUtils: we were expecting column name " + field + " but we saw " + strings[arrayIndex]);
                    } catch (IllegalArgumentException e) {
                        throw new RuntimeException("VCFReaderUtils: Unknown column name \"" + strings[arrayIndex] + "\", it does not match a known column header name.");
                    }
                    arrayIndex++;
                }
                while (arrayIndex < strings.length) {
                    if (!strings[arrayIndex].equals("FORMAT"))
                        auxTags.add(strings[arrayIndex]);
                    arrayIndex++;
                }
            } else {
                if ( str.startsWith("##INFO=") )
                    metaData.add(new VCFInfoHeaderLine(str.substring(7),version));
                else if ( str.startsWith("##FILTER=") )
                    metaData.add(new VCFFilterHeaderLine(str.substring(9),version));
                else if ( str.startsWith("##FORMAT=") )
                    metaData.add(new VCFFormatHeaderLine(str.substring(9),version));
                else {
                    int equals = str.indexOf("=");
                    if ( equals != -1 )
                        metaData.add(new VCFHeaderLine(str.substring(2, equals), str.substring(equals+1),version));
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
        return createRecord(line, mHeader, false);
    }

    public static VCFRecord createRecord(String line, VCFHeader mHeader, boolean ignoreGenotypes) {
        // things we need to make a VCF record
        Map<VCFHeader.HEADER_FIELDS, String> values = new HashMap<VCFHeader.HEADER_FIELDS, String>();
        String tokens[] = line.split("\\t");

        // check to ensure that the column count of tokens is right
        if (tokens.length != mHeader.getColumnCount()) {
            throw new RuntimeException("The input file line doesn't contain enough fields, it should have " + mHeader.getColumnCount() + " fields, it has " + tokens.length + ". Line = " + line);
        }

        int index = 0;
        for (VCFHeader.HEADER_FIELDS field : mHeader.getHeaderFields())
            values.put(field, tokens[index++]);
        // if we have genotyping data, we try and extract the genotype fields
        if ( ! ignoreGenotypes && mHeader.hasGenotypingData()) {
            String mFormatString = tokens[index];
            String keyStrings[] = mFormatString.split(":");
            List<VCFGenotypeRecord> genotypeRecords = new ArrayList<VCFGenotypeRecord>();
            index++;
			String[] alt_alleles = values.get(VCFHeader.HEADER_FIELDS.ALT).split(",");
            for (String str : mHeader.getGenotypeSamples()) {
                genotypeRecords.add(getVCFGenotype(str, keyStrings, tokens[index], alt_alleles, values.get(VCFHeader.HEADER_FIELDS.REF).charAt(0)));
                index++;
            }
            VCFRecord vrec = new VCFRecord(values, mFormatString, genotypeRecords);
            // associate the genotypes with this new record
            for ( VCFGenotypeRecord gr : genotypeRecords )
                gr.setVCFRecord(vrec);
            return vrec;

        }
        return new VCFRecord(values);
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
        return getVCFGenotype(sampleName, formatString.split(":"), genotypeString, altAlleles, referenceBase);
    }

    /**
     * generate a VCF genotype record, given it's format string, the genotype string, and allele info
     *
     * @param sampleName     the sample name
     * @param keyStrings     the split format string for this record, which contains the keys for the genotype parameters
     * @param genotypeString contains the phasing information, allele information, and values for genotype parameters
     * @param altAlleles     the alternate allele string array, which we index into based on the field parameters
     * @param referenceBase  the reference base
     *
     * @return a VCFGenotypeRecord
     */
    public static VCFGenotypeRecord getVCFGenotype(String sampleName, String[] keyStrings, String genotypeString, String altAlleles[], char referenceBase) {
        // parameters to create the VCF genotype record
        HashMap<String, String> tagToValue = new HashMap<String, String>();
        VCFGenotypeRecord.PHASE phase = VCFGenotypeRecord.PHASE.UNPHASED;
        List<VCFGenotypeEncoding> bases = new ArrayList<VCFGenotypeEncoding>();

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
            if (key.equals(VCFGenotypeRecord.GENOTYPE_KEY)) {
                Matcher m = gtPattern.matcher(parse);
                if (!m.matches())
                    throw new RuntimeException("VCFReaderUtils: Unable to match GT genotype flag to it's expected pattern, the field was: " + parse);
                phase = VCFGenotypeRecord.determinePhase(m.group(2));
                addAllele(m.group(1), altAlleles, referenceBase, bases);
                if (m.group(3).length() > 0) addAllele(m.group(3), altAlleles, referenceBase, bases);
            } else {
                if ( parse.length() == 0 )
                    parse = VCFGenotypeRecord.getMissingFieldValue(key);
                tagToValue.put(key, parse);
            }
            if (nextDivider + 1 >= genotypeString.length()) nextDivider = genotypeString.length() - 1;
            genotypeString = genotypeString.substring(nextDivider + 1, genotypeString.length());
        }
        if ( bases.size() > 0 && bases.get(0).equals(VCFGenotypeRecord.EMPTY_ALLELE) )
            tagToValue.clear();
        // catch some common errors, either there are too many field keys or there are two many field values
        else if ( keyStrings.length != tagToValue.size() + ((bases.size() > 0) ? 1 : 0))
            throw new RuntimeException("VCFReaderUtils: genotype value count doesn't match the key count (expected "
                    + keyStrings.length + " but saw " + tagToValue.size() + ")");
        else if ( genotypeString.length() > 0 )
            throw new RuntimeException("VCFReaderUtils: genotype string contained additional unprocessed fields: " + genotypeString
                    + ".  This most likely means that the format string is shorter then the value fields.");

        VCFGenotypeRecord rec = new VCFGenotypeRecord(sampleName, bases, phase);
        for ( Map.Entry<String, String> entry : tagToValue.entrySet() )
            rec.setField(entry.getKey(), entry.getValue());
        return rec;
    }


    /**
     * add an alternate allele to the list of alleles we have for a VCF genotype record
     *
     * @param alleleNumber  the allele number, as a string
     * @param altAlleles    the list of alternate alleles
     * @param referenceBase the reference base
     * @param bases         the list of bases for this genotype call
     */
    private static void addAllele(String alleleNumber, String[] altAlleles, char referenceBase, List<VCFGenotypeEncoding> bases) {
        if (alleleNumber.equals(VCFGenotypeRecord.EMPTY_ALLELE)) {
            bases.add(new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_ALLELE));
        } else {
            int alleleValue = Integer.valueOf(alleleNumber);
            // check to make sure the allele value is within bounds
            if (alleleValue < 0 || alleleValue > altAlleles.length)
                throw new IllegalArgumentException("VCFReaderUtils: the allele value of " + alleleValue + " is out of bounds given the alternate allele list.");
            if (alleleValue == 0)
                bases.add(new VCFGenotypeEncoding(String.valueOf(referenceBase)));
            else
                bases.add(new VCFGenotypeEncoding(altAlleles[alleleValue - 1]));
        }
    }
}
