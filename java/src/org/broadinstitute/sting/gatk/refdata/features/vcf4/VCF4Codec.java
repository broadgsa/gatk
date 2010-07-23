package org.broadinstitute.sting.gatk.refdata.features.vcf4;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.ParsingUtils;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.IOException;
import java.util.*;


/**
 * a feature codec for the VCF 4 specification.  Our aim is to read in the records and convert to VariantContext as
 * quickly as possible, relying on VariantContext to do the validation of any contradictory (or malformed) record parameters.
 */
public class VCF4Codec implements FeatureCodec, NameAwareCodec {


    // we have to store the list of strings that make up the header until they're needed
    private VCFHeader header = null;

    private VCFHeaderVersion version = VCFHeaderVersion.VCF4_0;
    // used to convert the index of the alternate allele in genotypes to a integer index
    private static int ZERO_CHAR = (byte)'0';

    // a mapping of the allele
    private static Map<String, List<Allele>> alleleMap = new HashMap<String, List<Allele>>(3);

    // cache the genotyope values
    private static String[] GTValueArray = new String[100];

    // for performance testing purposes
    public static boolean validate = true;

    // a key optimization -- we need a per thread string parts array, so we don't allocate a big array over and over
    // todo: make this thread safe?
    private String[] parts = null;

    // for performance we cache the hashmap of filter encodings for quick lookup
    private HashMap<String,LinkedHashSet<String>> filterHash = new HashMap<String,LinkedHashSet<String>>();

    // a set of the genotype keys?
    private String[] genotypeKeyArray = new String[100];

    // a mapping of the VCF fields to their type, filter fields, and format fields, for quick lookup to validate against
    TreeMap<String, VCFHeaderLineType> infoFields = new TreeMap<String, VCFHeaderLineType>();
    TreeMap<String, VCFHeaderLineType> formatFields = new TreeMap<String, VCFHeaderLineType>();
    ArrayList<String> filterFields = new ArrayList<String>();

    // do we want to validate the info, format, and filter fields
    private final boolean validateFromHeader = false;

    // we store a name to give to each of the variant contexts we emit
    private String name = "Unknown";

    private int lineNo = 0;

    // some classes need to transform the line before
    private LineTransform transformer = null;

    /**
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    @Override
    public Object readHeader(LineReader reader) {
        List<String> headerStrings = new ArrayList<String>();

        String line;
        try {
            boolean foundHeaderVersion = false;
            while ((line = reader.readLine()) != null) {
                lineNo++;
                if (line.startsWith(VCFHeader.METADATA_INDICATOR)) {
                    String[] lineFields = line.substring(2).split("=");
                    if (lineFields.length == 2 &&
                            VCFHeaderVersion.isVersionString(lineFields[1]) && VCFHeaderVersion.isFormatString(lineFields[0])) {
                        foundHeaderVersion = true;
                        this.version = VCFHeaderVersion.toHeaderVersion(lineFields[1]);
                    }
                    headerStrings.add(line);
                }
                else if (line.startsWith(VCFHeader.HEADER_INDICATOR)) {
                    if (!foundHeaderVersion) {
                        throw new CodecLineParsingException("We never saw a header line specifying VCF version");
                    }
                    return createHeader(headerStrings, line);
                }
                else {
                    throw new CodecLineParsingException("We never saw the required header line (starting with one #) for the input VCF file");
                }

            }
        } catch (IOException e) {
            throw new RuntimeException("IO Exception ", e);
        }
        throw new CodecLineParsingException("We never saw the required header line (starting with one #) for the input VCF file");

    }

    /**
     * create a VCF header
     * @param headerStrings a list of strings that represent all the ## entries
     * @param line the single # line (column names)
     * @return the count of header lines
     */
    public Object createHeader(List<String> headerStrings, String line) {
        headerStrings.add(line);
        header = VCFReaderUtils.createHeader(headerStrings, this.version);

        // setup our look-up lists for validation
        for ( VCFHeaderLine hl : header.getMetaData() ) {
            if ( hl instanceof VCFFilterHeaderLine )
                this.filterFields.add(((VCFFilterHeaderLine)hl).getName());
            if ( hl instanceof VCFFormatHeaderLine )
                this.formatFields.put(((VCFFormatHeaderLine)hl).getName(), ((VCFFormatHeaderLine)hl).getType());
            if ( hl instanceof VCFInfoHeaderLine )
                this.infoFields.put(((VCFInfoHeaderLine)hl).getName(), ((VCFInfoHeaderLine)hl).getType());
        }
        // sort the lists so we can binary search them later on
        Collections.sort(filterFields);

        return header;
    }

    /**
     * the fast decode function
     * @param line the line of text for the record
     * @return a feature, (not guaranteed complete) that has the correct start and stop
     */
    public Feature decodeLoc(String line) {
        return reallyDecode(line, false);
    }

    /**
     * decode the line into a feature (VariantContext)
     * @param line the line
     * @return a VariantContext
     */
    public Feature decode(String line) {
        return reallyDecode(line, true);
    }

    private Feature reallyDecode(String line, boolean parseGenotypes) {
        // the same line reader is not used for parsing the header and parsing lines, if we see a #, we've seen a header line
        if (line.startsWith(VCFHeader.HEADER_INDICATOR)) return null;

        if (parts == null)
            parts = new String[header.getColumnCount()];

        int nParts = ParsingUtils.split(line, parts, VCFConstants.FIELD_SEPARATOR.charAt(0));

        // our header cannot be null, we need the genotype sample names and counts
        if (header == null) throw new IllegalStateException("VCF Header cannot be null");

        // check to make sure the split resulted in the correct number of fields (8 + (1 + genotytpe counts if it has genotypes)
        if (nParts != header.getColumnCount())
            throw new IllegalArgumentException("we expected " + header.getColumnCount() + " columns and we got " + nParts + " for line " + line);

        return parseVCFLine(parts, parseGenotypes);
    }

    /**
     * create a an allele from an index and an array of alleles
     * @param index the index
     * @param alleles the alleles
     * @return an Allele
     */
    private static Allele oneAllele(char index, List<Allele> alleles) {
        if ( index == VCFConstants.EMPTY_ALLELE.charAt(0) )
            return Allele.NO_CALL;
        int i = ((byte)index) - ZERO_CHAR;
        return alleles.get(i);
    }


    /**
     * parse genotype alleles from the genotype string
     * @param GT         GT string
     * @param alleles    list of possible alleles
     * @param cache      cache of alleles for GT
     * @return the allele list for the GT string
     */
    private List<Allele> parseGenotypeAlleles(String GT, List<Allele> alleles, Map<String, List<Allele>> cache) {
        // this should cache results [since they are immutable] and return a single object for each genotype
        if ( GT.length() != 3 && GT.length() != 1 )
            throw new VCFParserException("Unreasonable number of alleles: " + "GT=" + GT + " length=" + GT.length()); // 0/1 => barf on 10/0

        List<Allele> GTAlleles = cache.get(GT);

        if ( GTAlleles == null ) {
            Allele allele1 = oneAllele(GT.charAt(0), alleles);
            GTAlleles = GT.length() == 3 ? Arrays.asList(allele1, oneAllele(GT.charAt(2), alleles)) : Arrays.asList(allele1);
            cache.put(GT, GTAlleles);
        }

        return GTAlleles;
    }

    /**
     * parse out the info fields
     * @param infoField the fields
     * @param id the indentifier
     * @return a mapping of keys to objects
     */
    private Map<String, Object> parseInfo(String infoField, String id) {
        Map<String, Object> attributes = new HashMap<String, Object>();

        if ( !infoField.equals(VCFConstants.EMPTY_INFO_FIELD) ) {
            for ( String field : Utils.split(infoField, VCFConstants.INFO_FIELD_SEPARATOR) ) {
                String key;
                Object value;

                int eqI = field.indexOf("=");
                if ( eqI != -1 ) {
                    key = field.substring(0, eqI);
                    String str = field.substring(eqI+1, field.length());

                    // lets see if the string contains a , separator
                    if ( str.contains(",") )
                        value = Arrays.asList(str.split(","));
                    else
                        value = str;
                } else {
                    key = field;
                    value = new Boolean(true);
                }

                attributes.put(key, value);
            }
        }
        // validate the fields
        validateFields(attributes.keySet(), new ArrayList<String>(infoFields.keySet()));

        attributes.put(VariantContext.ID_KEY, id);
        return attributes;
    }

    /**
     * validate the attributes against the stored fields of the appopriate type
     * @param attributes the list of fields to check for inclusion against the field array
     * @param fields the master list; all attributes must be in this list to validate
     */
    private void validateFields(Set<String> attributes, List<String> fields) {
        // validate the info fields
        if (validateFromHeader) {
            for (String attr : attributes)
                if (Collections.binarySearch(fields,attr) < 0)
                    throw new VCFParserException("Unable to find field describing attribute " + attr);
        }
    }

    /**
     * parse out the qual value
     * @param qualString the quality string
     * @return return a double
     */
    private Double parseQual(String qualString) {
        if ( qualString.equals(VCFConstants.MISSING_VALUE_v4) || qualString.equals(VCFConstants.MISSING_QUALITY_v3) )
            return VariantContext.NO_NEG_LOG_10PERROR;
        return Double.valueOf(qualString) / 10.0;
    }

    /**
     * parse out the alleles
     * @param ref the reference base
     * @param alts a string of alternates to break into alleles
     * @return a list of alleles, and a pair of the shortest and longest sequence
     */
    private List<Allele> parseAlleles(String ref, String alts) {
        List<Allele> alleles = new ArrayList<Allele>(2); // we are almost always biallelic
        // ref
        if (!checkAllele(ref, true))
            throw new VCFParserException("Unable to parse out correct reference allele, we saw = " + ref);
        Allele refAllele = Allele.create(ref, true);
        alleles.add(refAllele);

        if ( alts.indexOf(",") == -1 ) // only 1 alternatives, don't call string split
            parseSingleAllele(alleles, alts, false);
        else
            for ( String alt : Utils.split(alts, ",") )
                parseSingleAllele(alleles, alt, false);

        return alleles;
    }

    /**
     * check to make sure the allele is an acceptable allele
     * @param allele the allele to check
     * @param isRef are we the reference allele?
     * @return true if the allele is fine, false otherwise
     */
    private boolean checkAllele(String allele,boolean isRef) {
        if (allele.contains("<")) {
            Utils.warnUser("We are currently unable to parse out CNV encodings in VCF, we saw the following allele = " + allele);
            return false;
        }
        else if ( ! Allele.acceptableAlleleBases(allele,isRef) ) {
            throw new VCFParserException("Unparsable vcf record with allele " + allele);
        }
        return true;
    }

    /**
     * parse a single allele, given the allele list
     * @param alleles the alleles available
     * @param alt the allele to parse
     * @param isRef are we the reference allele?
     */
    private void parseSingleAllele(List<Allele> alleles, String alt, boolean isRef) {
        if (!checkAllele(alt,isRef))
            throw new VCFParserException("Unable to parse out correct alt allele, we saw = " + alt);

        Allele allele = Allele.create(alt, false);
        if ( ! allele.isNoCall() )
            alleles.add(allele);
    }

    /**
     * parse the filter string, first checking to see if we already have parsed it in a previous attempt
     * @param filterString the string to parse
     * @return a set of the filters applied
     */
    private Set<String> parseFilters(String filterString) {

        // null for unfiltered
        if ( filterString.equals(VCFConstants.UNFILTERED) )
            return null;

        // empty set for passes filters
        LinkedHashSet<String> fFields = new LinkedHashSet<String>();

        if ( this.version == VCFHeaderVersion.VCF4_0 ) {
            if ( filterString.equals(VCFConstants.PASSES_FILTERS_v4) )
                return fFields;
            if ( filterString.equals(VCFConstants.PASSES_FILTERS_v3) )
                throw new StingException(VCFConstants.PASSES_FILTERS_v3 + " is an invalid filter name in vcf4.0");
        } else if ( filterString.equals(VCFConstants.PASSES_FILTERS_v3) ) {
            return fFields;
        }

        // do we have the filter string cached?
        if ( filterHash.containsKey(filterString) )
            return filterHash.get(filterString);

        // otherwise we have to parse and cache the value
        if ( filterString.indexOf(VCFConstants.FILTER_CODE_SEPARATOR) == -1 )
            fFields.add(filterString);
        else
            fFields.addAll(Utils.split(filterString, VCFConstants.FILTER_CODE_SEPARATOR));

        filterHash.put(filterString, fFields);

        validateFields(fFields, filterFields);
        return fFields;
    }

    /**
     * parse out the VCF line
     *
     * @param parts the parts split up
     * @param parseGenotypes whether to parse genotypes or not
     * @return a variant context object
     */
    private VariantContext parseVCFLine(String[] parts, boolean parseGenotypes) {
//        try {
        // increment the line count
        lineNo++;

        // parse out the required fields
        String contig = parts[0];
        long pos = Long.valueOf(parts[1]);
        String id = parts[2];
        String ref = parts[3].toUpperCase();
        String alts = parts[4].toUpperCase();
        Double qual = parseQual(parts[5]);
        String filter = parts[6];
        String info = parts[7];

        // get our alleles, filters, and setup an attribute map
        List<Allele> alleles = parseAlleles(ref, alts);
        Set<String> filters = parseFilters(filter);
        Map<String, Object> attributes = parseInfo(info, id);

       // find out our current location, and clip the alleles down to their minimum length
        Pair<GenomeLoc, List<Allele>> locAndAlleles;
        if ( !isSingleNucleotideEvent(alleles) ) {
            if (this.version != VCFHeaderVersion.VCF4_0)
                throw new VCFParserException("Saw Indel/non SNP event on a VCF 3.3 or earlier file. Please convert file to VCF4.0 with VCFTools before using the GATK on it");
            locAndAlleles = clipAlleles(contig, pos, ref, alleles);
        } else {
            locAndAlleles = new Pair<GenomeLoc, List<Allele>>(GenomeLocParser.createGenomeLoc(contig, pos), alleles);
        }

        // a map to store our genotypes
        Map<String, Genotype> genotypes = null;

        // do we have genotyping data
        if (parts.length > 8 && parseGenotypes) {
            genotypes = createGenotypeMap(parts, locAndAlleles, 8);
        }

        VariantContext vc =  new VariantContext(name, locAndAlleles.first, locAndAlleles.second, genotypes, qual, filters, attributes);

        // Trim bases of all alleles if necessary
        return VariantContextUtils.createVariantContextWithTrimmedAlleles(vc);
    }

    private boolean isSingleNucleotideEvent(List<Allele> alleles) {
        for ( Allele a : alleles ) {
            if ( a.length() > 1 )
                return false;
        }
        return true;
    }

    class VCFParserException extends StingException {
        public VCFParserException(String msg) {
            super("Line " + lineNo + " generated parser exception " + msg);
        }

        public VCFParserException(String msg, Throwable throwable) {
            super("Line " + lineNo + " generated parser exception " + msg, throwable);
        }
    }

    /**
     * create a genotype map
     * @param parts the string parts
     * @param locAndAlleles the locations and the list of alleles
     * @param formatFieldLocation the position in the parts array that the genotype strings start
     * @return a mapping of sample name to genotype object
     */
    protected Map<String, Genotype> createGenotypeMap(String[] parts, Pair<GenomeLoc, List<Allele>> locAndAlleles, int formatFieldLocation) {
        Map<String, Genotype> genotypes = new LinkedHashMap<String, Genotype>(Math.max(parts.length - formatFieldLocation, 1));

        // get the format keys
        int nGTKeys = ParsingUtils.split(parts[formatFieldLocation], genotypeKeyArray, VCFConstants.GENOTYPE_FIELD_SEPARATOR.charAt(0));

        // cycle through the sample names
        Iterator<String> sampleNameIterator = header.getGenotypeSamples().iterator();

        // clear out our allele mapping
        alleleMap.clear();

        // cycle through the genotype strings
        for (int genotypeOffset = formatFieldLocation + 1; genotypeOffset < parts.length; genotypeOffset++) {
            int GTValueSplitSize = ParsingUtils.split(parts[genotypeOffset], GTValueArray, VCFConstants.GENOTYPE_FIELD_SEPARATOR.charAt(0));

            double GTQual = VariantContext.NO_NEG_LOG_10PERROR;
            Set<String> genotypeFilters = null;
            Map<String, String> gtAttributes = null;
            String sampleName = sampleNameIterator.next();

            // check to see if the value list is longer than the key list, which is a problem
            if (nGTKeys < GTValueSplitSize)
                throw new VCFParserException("Too few keys for compared to the value string " + sampleName + ", keys = " + parts[8] + " values = " + parts[genotypeOffset]);

            int genotypeAlleleLocation = -1;
            if (nGTKeys >= 1) {
                gtAttributes = new HashMap<String, String>(nGTKeys - 1);
                for (int i = 0; i < nGTKeys; i++) {
                    if (i >= GTValueSplitSize) {
                        if (genotypeKeyArray[i].equals(VCFConstants.GENOTYPE_QUALITY_KEY))
                            GTQual = parseQual(VCFConstants.MISSING_VALUE_v4);
                        else if (genotypeKeyArray[i].equals(VCFConstants.GENOTYPE_FILTER_KEY))
                            genotypeFilters = parseFilters(VCFConstants.MISSING_VALUE_v4);
                        else
                            gtAttributes.put(genotypeKeyArray[i],VCFConstants.MISSING_VALUE_v4);
                    }
                    else if (genotypeKeyArray[i].equals(VCFConstants.GENOTYPE_KEY))
                        if (i != 0)
                            throw new VCFParserException("Saw GT at position " + i + ", it must be at the first position for genotypes. At location = " + locAndAlleles.first);
                        else
                            genotypeAlleleLocation = i;
                    else if (genotypeKeyArray[i].equals(VCFConstants.GENOTYPE_QUALITY_KEY))
                        GTQual = parseQual(GTValueArray[i]);
                    else if (genotypeKeyArray[i].equals(VCFConstants.GENOTYPE_FILTER_KEY))
                        genotypeFilters = parseFilters(GTValueArray[i]);
                    else {
                        if (this.version != VCFHeaderVersion.VCF4_0 && GTValueArray[i].equals(VCFConstants.MISSING_GENOTYPE_QUALITY_v3))
                            GTValueArray[i] = VCFConstants.MISSING_VALUE_v4;
                        gtAttributes.put(genotypeKeyArray[i], GTValueArray[i]);
                    }
                }
                // validate the format fields
                validateFields(gtAttributes.keySet(), new ArrayList<String>(formatFields.keySet()));
            }
            // check to make sure we found a gentoype field
            if (genotypeAlleleLocation < 0) throw new VCFParserException("Unable to find required field GT for record " + locAndAlleles.first);

            // assuming allele list length in the single digits, could be bad.  Check for > 1 for haploid genotypes
            boolean phased = GTValueArray[genotypeAlleleLocation].length() > 1 && GTValueArray[genotypeAlleleLocation].charAt(1) == '|';

            // add it to the list
            genotypes.put(sampleName, new Genotype(sampleName,
                    parseGenotypeAlleles(GTValueArray[genotypeAlleleLocation], locAndAlleles.second, alleleMap),
                    GTQual,
                    genotypeFilters,
                    gtAttributes,
                    phased));

        }
        return genotypes;
    }

    /**
     * clip the alleles, based on the reference
     *
     * @param contig our contig position
     * @param position the unadjusted start position (pre-clipping)
     * @param ref the reference string
     * @param unclippedAlleles the list of unclipped alleles
     * @return a list of alleles, clipped to the reference
     */
    static Pair<GenomeLoc,List<Allele>> clipAlleles(String contig, long position, String ref, List<Allele> unclippedAlleles) {
        List<Allele> newAlleleList = new ArrayList<Allele>();

        // find the preceeding string common to all alleles and the reference
        boolean clipping = true;
        for (Allele a : unclippedAlleles)
                if (a.length() < 1 || (a.getBases()[0] != ref.getBytes()[0])) {
                    clipping = false;
                }
        int forwardClipping = (clipping) ? 1 : 0;

        int reverseClipped = 0;
        clipping = true;
        while (clipping) {
            for (Allele a : unclippedAlleles)
                if (a.length() - reverseClipped <= forwardClipping || a.length() - forwardClipping == 0)
                    clipping = false;
                else if (a.getBases()[a.length()-reverseClipped-1] != ref.getBytes()[ref.length()-reverseClipped-1])
                    clipping = false;
            if (clipping) reverseClipped++;
        }

        for (Allele a : unclippedAlleles)
            newAlleleList.add(Allele.create(Arrays.copyOfRange(a.getBases(),forwardClipping,a.getBases().length-reverseClipped),a.isReference()));

        // the new reference length
        int refLength = ref.length() - reverseClipped;

        return new Pair<GenomeLoc,List<Allele>>(GenomeLocParser.createGenomeLoc(contig,position,(position+Math.max(refLength - 1,0))),
                newAlleleList);
    }

    /**
     *
     * @return the type of record
     */
    @Override
    public Class getFeatureType() {
        return VariantContext.class;
    }

    /**
     * get the name of this codec
     * @return our set name
     */
    public String getName() {
        return name;
    }

    /**
     * set the name of this codec
     * @param name new name
     */
    public void setName(String name) {
        this.name = name;
    }

    public static interface LineTransform {
        public String lineTransform(String line);
    }

    public LineTransform getTransformer() {
        return transformer;
    }

    public void setTransformer(LineTransform transformer) {
        this.transformer = transformer;
    }

}
