package org.broadinstitute.sting.gatk.refdata.features.vcf4;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.util.LineReader;
import org.broad.tribble.util.ParsingUtils;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
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

    // a variant context flag for original allele strings
    public static final String ORIGINAL_ALLELE_LIST = "ORIGINAL_ALLELE_LIST";

    // we have to store the list of strings that make up the header until they're needed
    private VCFHeader header = null;

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

    /**
     * this method is a big hack, since I haven't gotten to updating the VCF header for the 4.0 updates
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    @Override
    public int readHeader(LineReader reader) {
        List<String> headerStrings = new ArrayList<String>();

        String line = "";
        try {
            while ((line = reader.readLine()) != null) {
                lineNo++;
                if (line.startsWith("##")) {
                    headerStrings.add(line);
                }
                else if (line.startsWith("#")) {
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
    public int createHeader(List<String> headerStrings, String line) {
        headerStrings.add(line);
        header = VCFReaderUtils.createHeader(headerStrings, VCFHeaderVersion.VCF4_0);

        // load the parsing fields
        Set<VCFHeaderLine> headerLines = header.getMetaData();

        // setup our look-up lists for validation
        for (VCFHeaderLine hl : headerLines) {
            if (hl.getClass() == VCFFilterHeaderLine.class)
                this.filterFields.add(((VCFFilterHeaderLine)hl).getName());
            if (hl.getClass() == VCFFormatHeaderLine.class)
                                        this.formatFields.put(((VCFFormatHeaderLine)hl).getName(),((VCFFormatHeaderLine)hl).getType());
            if (hl.getClass() == VCFInfoHeaderLine.class)
                                        this.infoFields.put(((VCFInfoHeaderLine)hl).getName(),((VCFInfoHeaderLine)hl).getType());
        }
        // sort the lists so we can binary search them later on
        Collections.sort(filterFields);

        return headerStrings.size();
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
        if (parts == null)
            parts = new String[header.getColumnCount()];

        int nParts = ParsingUtils.split(line, parts, '\t');

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
        if ( index == '.' )
            return Allele.NO_CALL;
        else {
            int i = ((byte)index) - ZERO_CHAR;
            return alleles.get(i);
        }
    }


    /**
     * parse genotype alleles from the genotype string
     * @param GT
     * @param alleles
     * @param cache
     * @return
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

        if ( ! infoField.equals(".") ) { // empty info field
            for ( String field : Utils.split(infoField, ";") ) {
                String key;
                Object value;

                int eqI = field.indexOf("=");
                if ( eqI != -1 ) {
                    key = field.substring(0, eqI);
                    String str = field.substring(eqI+1, field.length());

                    // lets see if the string contains a , separator
                    if (str.contains(",")) {
                        List<Object> objects = new ArrayList<Object>();
                        String[] split = str.split(",");
                        for (String substring : split) {
                            VCFHeaderLineType type = infoFields.get(key);
                            objects.add(type != null ? type.convert(substring,VCFCompoundHeaderLine.SupportedHeaderLineType.INFO) : substring);
                        }
                        value = objects;
                    } else {
                        VCFHeaderLineType type = infoFields.get(key);
                        value = type != null ? type.convert(str,VCFCompoundHeaderLine.SupportedHeaderLineType.INFO) : str;
                    }
                    //System.out.printf("%s %s%n", key, value);
                } else {
                    key = field;
                    value = 1;
                }

                attributes.put(key, value);
            }
        }
        // validate the fields
        validateFields(attributes.keySet(),new ArrayList(infoFields.keySet()));

        attributes.put("ID", id);
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
            int count = 0;
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
    private static Double parseQual(String qualString) {
        // todo -- remove double once we deal with annoying VCFs from 1KG
        return qualString.equals("-1") || qualString.equals("-1.0") ? VariantContext.NO_NEG_LOG_10PERROR : Double.valueOf(qualString) / 10;
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
        Set<String> fFields;

        // a PASS is simple (no filters)
        if ( filterString.equals("PASS") ) {
            return null;
        }
        // else do we have the filter string cached?
        else if (filterHash.containsKey(filterString)) {
            fFields = filterHash.get(filterString);
        }
        // otherwise we have to parse and cache the value
        else {
            LinkedHashSet<String> s = new LinkedHashSet<String>(1);
            if ( filterString.indexOf(";") == -1 ) {
                s.add(filterString);
            } else {
                s.addAll(Utils.split(filterString, ";"));
            }
            filterHash.put(filterString,s);
            fFields = s;
        }

        validateFields(fFields,filterFields);
        return fFields;
    }

    /**
     * parse out the VCF line
     *
     * @param parts the parts split up
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
            if ( hasIndel(alleles) ) {
                attributes.put(ORIGINAL_ALLELE_LIST,alleles);
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

            return new VariantContext(name, locAndAlleles.first, locAndAlleles.second, genotypes, qual, filters, attributes);
    }

    private boolean hasIndel(List<Allele> alleles) {
        int lengthOfFirstEntry = alleles.get(0).length();
        for ( int i = 1; i < alleles.size(); i++ ) {
            if ( alleles.get(i).length() != lengthOfFirstEntry )
                return true;
        }
        return false;
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
        int nGTKeys = ParsingUtils.split(parts[formatFieldLocation], genotypeKeyArray, ':');

        // cycle through the sample names
        Iterator<String> sampleNameIterator = header.getGenotypeSamples().iterator();

        // clear out our allele mapping
        alleleMap.clear();

        // cycle through the genotype strings
        for (int genotypeOffset = formatFieldLocation + 1; genotypeOffset < parts.length; genotypeOffset++) {
            int GTValueSplitSize = ParsingUtils.split(parts[genotypeOffset], GTValueArray, ':');

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
                    if (i >= GTValueSplitSize)
                        gtAttributes.put(genotypeKeyArray[i],".");
                    else if (genotypeKeyArray[i].equals("GT"))
                        if (i != 0)
                            throw new VCFParserException("Saw GT at position " + i + ", it must be at the first position for genotypes. At location = " + locAndAlleles.first);
                        else
                            genotypeAlleleLocation = i;
                    else if (genotypeKeyArray[i].equals("GQ"))
                        GTQual = parseQual(GTValueArray[i]);
                    else if (genotypeKeyArray[i].equals("FT")) // deal with genotype filters here
                        genotypeFilters = parseFilters(GTValueArray[i]);
                    else
                        gtAttributes.put(genotypeKeyArray[i], GTValueArray[i]);

                }
                // validate the format fields
                validateFields(gtAttributes.keySet(), new ArrayList(formatFields.keySet()));
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
        int refLength = ref.length() - forwardClipping - reverseClipped;

        return new Pair<GenomeLoc,List<Allele>>(GenomeLocParser.createGenomeLoc(contig,position+forwardClipping,(position+forwardClipping+Math.max(refLength - 1,0))),
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
     * get the header
     * @param clazz the class were expecting
     * @return our VCFHeader
     * @throws ClassCastException
     */
    @Override
    public VCFHeader getHeader(Class clazz) throws ClassCastException {
        if (clazz != VCFHeader.class) throw new ClassCastException("expecting class " + clazz + " but VCF4Codec provides " + VCFHeader.class);
        return this.header;
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
     * @param name
     */
    public void setName(String name) {
        this.name = name;
    }
}
