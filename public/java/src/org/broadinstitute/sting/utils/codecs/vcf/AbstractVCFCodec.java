package org.broadinstitute.sting.utils.codecs.vcf;

import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.NameAwareCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.BlockCompressedInputStream;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.gatk.refdata.SelfScopingFeatureCodec;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;


public abstract class AbstractVCFCodec implements FeatureCodec, NameAwareCodec, VCFParser, SelfScopingFeatureCodec {

    protected final static Logger log = Logger.getLogger(VCFCodec.class);
    protected final static int NUM_STANDARD_FIELDS = 8;  // INFO is the 8th column

    protected VCFHeaderVersion version;

    // we have to store the list of strings that make up the header until they're needed
    protected VCFHeader header = null;

    // a mapping of the allele
    protected Map<String, List<Allele>> alleleMap = new HashMap<String, List<Allele>>(3);

    // for ParsingUtils.split
    protected String[] GTValueArray = new String[100];
    protected String[] genotypeKeyArray = new String[100];
    protected String[] infoFieldArray = new String[1000];
    protected String[] infoValueArray = new String[1000];

    // for performance testing purposes
    public static boolean validate = true;

    // a key optimization -- we need a per thread string parts array, so we don't allocate a big array over and over
    // todo: make this thread safe?
    protected String[] parts = null;
    protected String[] genotypeParts = null;

    // for performance we cache the hashmap of filter encodings for quick lookup
    protected HashMap<String,LinkedHashSet<String>> filterHash = new HashMap<String,LinkedHashSet<String>>();

    // a mapping of the VCF fields to their type, filter fields, and format fields, for quick lookup to validate against
    TreeMap<String, VCFHeaderLineType> infoFields = new TreeMap<String, VCFHeaderLineType>();
    TreeMap<String, VCFHeaderLineType> formatFields = new TreeMap<String, VCFHeaderLineType>();
    Set<String> filterFields = new HashSet<String>();

    // we store a name to give to each of the variant contexts we emit
    protected String name = "Unknown";

    protected int lineNo = 0;

    protected Map<String, String> stringCache = new HashMap<String, String>();


    /**
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    public abstract Object readHeader(LineReader reader);

    /**
     * create a genotype map
     * @param str the string
     * @param alleles the list of alleles
     * @param chr chrom
     * @param pos position
     * @return a mapping of sample name to genotype object
     */
    public abstract Map<String, Genotype> createGenotypeMap(String str, List<Allele> alleles, String chr, int pos);


    /**
     * parse the filter string, first checking to see if we already have parsed it in a previous attempt
     * @param filterString the string to parse
     * @return a set of the filters applied
     */
    protected abstract Set<String> parseFilters(String filterString);

    /**
     * create a VCF header
     * @param headerStrings a list of strings that represent all the ## entries
     * @param line the single # line (column names)
     * @return the count of header lines
     */
    protected Object createHeader(List<String> headerStrings, String line) {

        headerStrings.add(line);

        Set<VCFHeaderLine> metaData = new TreeSet<VCFHeaderLine>();
        Set<String> auxTags = new LinkedHashSet<String>();
        // iterate over all the passed in strings
        for ( String str : headerStrings ) {
            if ( !str.startsWith(VCFHeader.METADATA_INDICATOR) ) {
                String[] strings = str.substring(1).split(VCFConstants.FIELD_SEPARATOR);
                if ( strings.length < VCFHeader.HEADER_FIELDS.values().length )
                    throw new TribbleException.InvalidHeader("there are not enough columns present in the header line: " + str);

                int arrayIndex = 0;
                for (VCFHeader.HEADER_FIELDS field : VCFHeader.HEADER_FIELDS.values()) {
                    try {
                        if (field != VCFHeader.HEADER_FIELDS.valueOf(strings[arrayIndex]))
                            throw new TribbleException.InvalidHeader("we were expecting column name '" + field + "' but we saw '" + strings[arrayIndex] + "'");
                    } catch (IllegalArgumentException e) {
                        throw new TribbleException.InvalidHeader("unknown column name '" + strings[arrayIndex] + "'; it does not match a legal column header name.");
                    }
                    arrayIndex++;
                }

                boolean sawFormatTag = false;
                if ( arrayIndex < strings.length ) {
                    if ( !strings[arrayIndex].equals("FORMAT") )
                        throw new TribbleException.InvalidHeader("we were expecting column name 'FORMAT' but we saw '" + strings[arrayIndex] + "'");
                    sawFormatTag = true;
                    arrayIndex++;
                }

                while ( arrayIndex < strings.length )
                    auxTags.add(strings[arrayIndex++]);

                if ( sawFormatTag && auxTags.size() == 0 )
                    throw new UserException.MalformedVCFHeader("The FORMAT field was provided but there is no genotype/sample data");

            } else {
                if ( str.startsWith("##INFO=") ) {
                    VCFInfoHeaderLine info = new VCFInfoHeaderLine(str.substring(7),version);
                    metaData.add(info);
                    infoFields.put(info.getName(), info.getType());
                } else if ( str.startsWith("##FILTER=") ) {
                    VCFFilterHeaderLine filter = new VCFFilterHeaderLine(str.substring(9),version);
                    metaData.add(filter);
                    filterFields.add(filter.getName());
                } else if ( str.startsWith("##FORMAT=") ) {
                    VCFFormatHeaderLine format = new VCFFormatHeaderLine(str.substring(9),version);
                    metaData.add(format);
                    formatFields.put(format.getName(), format.getType());
                } else {
                    int equals = str.indexOf("=");
                    if ( equals != -1 )
                        metaData.add(new VCFHeaderLine(str.substring(2, equals), str.substring(equals+1)));
                }
            }
        }

        header = new VCFHeader(metaData, auxTags);
        return header;
    }

    /**
     * the fast decode function
     * @param line the line of text for the record
     * @return a feature, (not guaranteed complete) that has the correct start and stop
     */
    public Feature decodeLoc(String line) {
        lineNo++;

        // the same line reader is not used for parsing the header and parsing lines, if we see a #, we've seen a header line
        if (line.startsWith(VCFHeader.HEADER_INDICATOR)) return null;

        // our header cannot be null, we need the genotype sample names and counts
        if (header == null) throw new ReviewedStingException("VCF Header cannot be null when decoding a record");

        final String[] locParts = new String[6];
        int nParts = ParsingUtils.split(line, locParts, VCFConstants.FIELD_SEPARATOR_CHAR, true);

        if ( nParts != 6 )
            throw new UserException.MalformedVCF("there aren't enough columns for line " + line, lineNo);

        // get our alleles (because the end position depends on them)
        final String ref = getCachedString(locParts[3].toUpperCase());
        final String alts = getCachedString(locParts[4].toUpperCase());
        final List<Allele> alleles = parseAlleles(ref, alts, lineNo);

        // find out our location
        final int start = Integer.valueOf(locParts[1]);
        int stop = start;

        // ref alleles don't need to be single bases for monomorphic sites
        if ( alleles.size() == 1 ) {
            stop = start + alleles.get(0).length() - 1;
        } else if ( !isSingleNucleotideEvent(alleles) ) {
            stop = clipAlleles(start, ref, alleles, null, lineNo);
        }

        return new VCFLocFeature(locParts[0], start, stop);
    }

    private final static class VCFLocFeature implements Feature {

        final String chr;
        final int start, stop;

        private VCFLocFeature(String chr, int start, int stop) {
            this.chr = chr;
            this.start = start;
            this.stop = stop;
        }

        public String getChr() { return chr; }
        public int getStart() { return start; }
        public int getEnd() { return stop; }
    }


    /**
     * decode the line into a feature (VariantContext)
     * @param line the line
     * @return a VariantContext
     */
    public Feature decode(String line) {
        // the same line reader is not used for parsing the header and parsing lines, if we see a #, we've seen a header line
        if (line.startsWith(VCFHeader.HEADER_INDICATOR)) return null;

        // our header cannot be null, we need the genotype sample names and counts
        if (header == null) throw new ReviewedStingException("VCF Header cannot be null when decoding a record");

        if (parts == null)
            parts = new String[Math.min(header.getColumnCount(), NUM_STANDARD_FIELDS+1)];

        int nParts = ParsingUtils.split(line, parts, VCFConstants.FIELD_SEPARATOR_CHAR, true);

        // if we have don't have a header, or we have a header with no genotyping data check that we have eight columns.  Otherwise check that we have nine (normal colummns + genotyping data)
        if (( (header == null || !header.hasGenotypingData()) && nParts != NUM_STANDARD_FIELDS) ||
             (header != null && header.hasGenotypingData() && nParts != (NUM_STANDARD_FIELDS + 1)) )
            throw new UserException.MalformedVCF("there aren't enough columns for line " + line + " (we expected " + (header == null ? NUM_STANDARD_FIELDS : NUM_STANDARD_FIELDS + 1) +
                    " tokens, and saw " + nParts + " )", lineNo);

        return parseVCFLine(parts);
    }

    protected void generateException(String message) {
        throw new UserException.MalformedVCF(message, lineNo);
    }

    protected static void generateException(String message, int lineNo) {
        throw new UserException.MalformedVCF(message, lineNo);
    }

    /**
     * parse out the VCF line
     *
     * @param parts the parts split up
     * @return a variant context object
     */
    private VariantContext parseVCFLine(String[] parts) {
        // increment the line count
        lineNo++;

        // parse out the required fields
        String contig = getCachedString(parts[0]);
        int pos = Integer.valueOf(parts[1]);
        String id = null;
        if ( parts[2].length() == 0 )
            generateException("The VCF specification requires a valid ID field");
        else if ( parts[2].equals(VCFConstants.EMPTY_ID_FIELD) )
            id = VCFConstants.EMPTY_ID_FIELD;
        else
            id = new String(parts[2]);
        String ref = getCachedString(parts[3].toUpperCase());
        String alts = getCachedString(parts[4].toUpperCase());
        Double qual = parseQual(parts[5]);
        String filter = getCachedString(parts[6]);
        String info = new String(parts[7]);

        // get our alleles, filters, and setup an attribute map
        List<Allele> alleles = parseAlleles(ref, alts, lineNo);
        Set<String> filters = parseFilters(filter);
        Map<String, Object> attributes = parseInfo(info, id);

        // find out our current location, and clip the alleles down to their minimum length
        int loc = pos;
        // ref alleles don't need to be single bases for monomorphic sites
        if ( alleles.size() == 1 ) {
            loc = pos + alleles.get(0).length() - 1;
        } else if ( !isSingleNucleotideEvent(alleles) ) {
            ArrayList<Allele> newAlleles = new ArrayList<Allele>();
            loc = clipAlleles(pos, ref, alleles, newAlleles, lineNo);
            alleles = newAlleles;
        }

        // do we have genotyping data
        if (parts.length > NUM_STANDARD_FIELDS) {
            attributes.put(VariantContext.UNPARSED_GENOTYPE_MAP_KEY, new String(parts[8]));
            attributes.put(VariantContext.UNPARSED_GENOTYPE_PARSER_KEY, this);
        }

        VariantContext vc = null;
        try {
            vc =  new VariantContext(name, contig, pos, loc, alleles, qual, filters, attributes, ref.getBytes()[0]);
        } catch (Exception e) {
            generateException(e.getMessage());
        }

        // did we resort the sample names?  If so, we need to load the genotype data
        if ( !header.samplesWereAlreadySorted() )
            vc.getGenotypes();

        return vc;
    }

    /**
     *
     * @return the type of record
     */
    public Class<VariantContext> getFeatureType() {
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

    /**
     * Return a cached copy of the supplied string.
     *
     * @param str string
     * @return interned string
     */
    protected String getCachedString(String str) {
        String internedString = stringCache.get(str);
        if ( internedString == null ) {
            internedString = new String(str);
            stringCache.put(internedString, internedString);
        }
        return internedString;
    }

    /**
     * parse out the info fields
     * @param infoField the fields
     * @param id the indentifier
     * @return a mapping of keys to objects
     */
    private Map<String, Object> parseInfo(String infoField, String id) {
        Map<String, Object> attributes = new HashMap<String, Object>();

        if ( infoField.length() == 0 )
            generateException("The VCF specification requires a valid info field");

        if ( !infoField.equals(VCFConstants.EMPTY_INFO_FIELD) ) {
            if ( infoField.indexOf("\t") != -1 || infoField.indexOf(" ") != -1 )
                generateException("The VCF specification does not allow for whitespace in the INFO field");

            int infoFieldSplitSize = ParsingUtils.split(infoField, infoFieldArray, VCFConstants.INFO_FIELD_SEPARATOR_CHAR, false);
            for (int i = 0; i < infoFieldSplitSize; i++) {
                String key;
                Object value;

                int eqI = infoFieldArray[i].indexOf("=");
                if ( eqI != -1 ) {
                    key = infoFieldArray[i].substring(0, eqI);
                    String str = infoFieldArray[i].substring(eqI+1);

                    // split on the INFO field separator
                    int infoValueSplitSize = ParsingUtils.split(str, infoValueArray, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR_CHAR, false);
                    if ( infoValueSplitSize == 1 ) {
                        value = infoValueArray[0];
                    } else {
                        ArrayList<String> valueList = new ArrayList<String>(infoValueSplitSize);
                        for ( int j = 0; j < infoValueSplitSize; j++ )
                            valueList.add(infoValueArray[j]);
                        value = valueList;
                    }
                } else {
                    key = infoFieldArray[i];
                    value = true;
                }

                attributes.put(key, value);
            }
        }

        if ( ! id.equals(VCFConstants.EMPTY_ID_FIELD) )
            attributes.put(VariantContext.ID_KEY, id);
        return attributes;
    }

    /**
     * create a an allele from an index and an array of alleles
     * @param index the index
     * @param alleles the alleles
     * @return an Allele
     */
    protected static Allele oneAllele(String index, List<Allele> alleles) {
        if ( index.equals(VCFConstants.EMPTY_ALLELE) )
            return Allele.NO_CALL;
        int i = Integer.valueOf(index);
        if ( i >= alleles.size() )
            throw new TribbleException.InternalCodecException("The allele with index " + index + " is not defined in the REF/ALT columns in the record");
        return alleles.get(i);
    }


    /**
     * parse genotype alleles from the genotype string
     * @param GT         GT string
     * @param alleles    list of possible alleles
     * @param cache      cache of alleles for GT
     * @return the allele list for the GT string
     */
    protected static List<Allele> parseGenotypeAlleles(String GT, List<Allele> alleles, Map<String, List<Allele>> cache) {
        // cache results [since they are immutable] and return a single object for each genotype
        List<Allele> GTAlleles = cache.get(GT);

        if ( GTAlleles == null ) {
            StringTokenizer st = new StringTokenizer(GT, VCFConstants.PHASING_TOKENS);
            GTAlleles = new ArrayList<Allele>(st.countTokens());
            while ( st.hasMoreTokens() ) {
                String genotype = st.nextToken();
                GTAlleles.add(oneAllele(genotype, alleles));
            }
            cache.put(GT, GTAlleles);
        }

        return GTAlleles;
    }

    /**
     * parse out the qual value
     * @param qualString the quality string
     * @return return a double
     */
    protected static Double parseQual(String qualString) {
        // if we're the VCF 4 missing char, return immediately
        if ( qualString.equals(VCFConstants.MISSING_VALUE_v4))
            return VariantContext.NO_NEG_LOG_10PERROR;

        Double val = Double.valueOf(qualString);

        // check to see if they encoded the missing qual score in VCF 3 style, with either the -1 or -1.0.  check for val < 0 to save some CPU cycles
        if ((val < 0) && (Math.abs(val - VCFConstants.MISSING_QUALITY_v3_DOUBLE) < VCFConstants.VCF_ENCODING_EPSILON))
            return VariantContext.NO_NEG_LOG_10PERROR;

        // scale and return the value
        return val / 10.0;
    }

    /**
     * parse out the alleles
     * @param ref the reference base
     * @param alts a string of alternates to break into alleles
     * @param lineNo  the line number for this record
     * @return a list of alleles, and a pair of the shortest and longest sequence
     */
    protected static List<Allele> parseAlleles(String ref, String alts, int lineNo) {
        List<Allele> alleles = new ArrayList<Allele>(2); // we are almost always biallelic
        // ref
        checkAllele(ref, true, lineNo);
        Allele refAllele = Allele.create(ref, true);
        alleles.add(refAllele);

        if ( alts.indexOf(",") == -1 ) // only 1 alternatives, don't call string split
            parseSingleAltAllele(alleles, alts, lineNo);
        else
            for ( String alt : alts.split(",") )
                parseSingleAltAllele(alleles, alt, lineNo);

        return alleles;
    }

    /**
     * check to make sure the allele is an acceptable allele
     * @param allele the allele to check
     * @param isRef are we the reference allele?
     * @param lineNo  the line number for this record
     */
    private static void checkAllele(String allele, boolean isRef, int lineNo) {
	if ( allele == null || allele.length() == 0 )
	    generateException("Empty alleles are not permitted in VCF records", lineNo);

        if ( isSymbolicAllele(allele) ) {
            if ( isRef ) {
                generateException("Symbolic alleles not allowed as reference allele: " + allele, lineNo);
            }
        } else {
            // check for VCF3 insertions or deletions
            if ( (allele.charAt(0) == VCFConstants.DELETION_ALLELE_v3) || (allele.charAt(0) == VCFConstants.INSERTION_ALLELE_v3) )
                generateException("Insertions/Deletions are not supported when reading 3.x VCF's. Please" +
                        " convert your file to VCF4 using VCFTools, available at http://vcftools.sourceforge.net/index.html", lineNo);

            if (!Allele.acceptableAlleleBases(allele))
                generateException("Unparsable vcf record with allele " + allele, lineNo);

            if ( isRef && allele.equals(VCFConstants.EMPTY_ALLELE) )
                generateException("The reference allele cannot be missing", lineNo);
        }
    }

    /**
     * return true if this is a symbolic allele (e.g. <SOMETAG>) otherwise false
     * @param allele the allele to check
     * @return true if the allele is a symbolic allele, otherwise false
     */
    private static boolean isSymbolicAllele(String allele) {
        return (allele != null && allele.startsWith("<") && allele.endsWith(">") && allele.length() > 2);
    }

    /**
     * parse a single allele, given the allele list
     * @param alleles the alleles available
     * @param alt the allele to parse
     * @param lineNo  the line number for this record
     */
    private static void parseSingleAltAllele(List<Allele> alleles, String alt, int lineNo) {
        checkAllele(alt, false, lineNo);

        Allele allele = Allele.create(alt, false);
        if ( ! allele.isNoCall() )
            alleles.add(allele);
    }

    protected static boolean isSingleNucleotideEvent(List<Allele> alleles) {
        for ( Allele a : alleles ) {
            if ( a.length() != 1 )
                return false;
        }
        return true;
    }

    public static int computeForwardClipping(List<Allele> unclippedAlleles, String ref) {
        boolean clipping = true;

        for ( Allele a : unclippedAlleles ) {
            if ( a.isSymbolic() )
                continue;

            if ( a.length() < 1 || (a.getBases()[0] != ref.getBytes()[0]) ) {
                clipping = false;
                break;
            }
        }

        return (clipping) ? 1 : 0;
    }

    protected static int computeReverseClipping(List<Allele> unclippedAlleles, String ref, int forwardClipping, int lineNo) {
        int clipping = 0;
        boolean stillClipping = true;

        while ( stillClipping ) {
            for ( Allele a : unclippedAlleles ) {
                if ( a.isSymbolic() )
                    continue;

                if ( a.length() - clipping <= forwardClipping || a.length() - forwardClipping == 0 )
                    stillClipping = false;
                else if ( ref.length() == clipping )
                    generateException("bad alleles encountered", lineNo);
                else if ( a.getBases()[a.length()-clipping-1] != ref.getBytes()[ref.length()-clipping-1] )
                    stillClipping = false;
            }
            if ( stillClipping )
                clipping++;
        }

        return clipping;
    }
    /**
     * clip the alleles, based on the reference
     *
     * @param position the unadjusted start position (pre-clipping)
     * @param ref the reference string
     * @param unclippedAlleles the list of unclipped alleles
     * @param clippedAlleles output list of clipped alleles
     * @param lineNo the current line number in the file
     * @return the new reference end position of this event
     */
    protected static int clipAlleles(int position, String ref, List<Allele> unclippedAlleles, List<Allele> clippedAlleles, int lineNo) {

        int forwardClipping = computeForwardClipping(unclippedAlleles, ref);
        int reverseClipping = computeReverseClipping(unclippedAlleles, ref, forwardClipping, lineNo);

        if ( clippedAlleles != null ) {
            for ( Allele a : unclippedAlleles ) {
                if ( a.isSymbolic() ) {
                    clippedAlleles.add(a);
                } else {
                    clippedAlleles.add(Allele.create(Arrays.copyOfRange(a.getBases(), forwardClipping, a.getBases().length-reverseClipping), a.isReference()));
                }
            }
        }

        // the new reference length
        int refLength = ref.length() - reverseClipping;

        return position+Math.max(refLength - 1,0);
    }

    public final static boolean canDecodeFile(final File potentialInput, final String MAGIC_HEADER_LINE) {
        try {
            return isVCFStream(new FileInputStream(potentialInput), MAGIC_HEADER_LINE) ||
                    isVCFStream(new GZIPInputStream(new FileInputStream(potentialInput)), MAGIC_HEADER_LINE) ||
                    isVCFStream(new BlockCompressedInputStream(new FileInputStream(potentialInput)), MAGIC_HEADER_LINE);
        } catch ( FileNotFoundException e ) {
            return false;
        } catch ( IOException e ) {
            return false;
        }
    }

    private final static boolean isVCFStream(final InputStream stream, final String MAGIC_HEADER_LINE) {
        try {
            byte[] buff = new byte[MAGIC_HEADER_LINE.length()];
            int nread = stream.read(buff, 0, MAGIC_HEADER_LINE.length());
            boolean eq = Arrays.equals(buff, MAGIC_HEADER_LINE.getBytes());
            return eq;
//            String firstLine = new String(buff);
//            return firstLine.startsWith(MAGIC_HEADER_LINE);
        } catch ( IOException e ) {
            return false;
        } catch ( RuntimeException e ) {
            return false;
        } finally {
            try { stream.close(); } catch ( IOException e ) {}
        }
    }
}
