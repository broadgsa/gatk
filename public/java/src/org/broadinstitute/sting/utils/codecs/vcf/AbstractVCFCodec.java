package org.broadinstitute.sting.utils.codecs.vcf;

import org.apache.log4j.Logger;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.NameAwareCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.BlockCompressedInputStream;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.zip.GZIPInputStream;


public abstract class AbstractVCFCodec extends AsciiFeatureCodec<VariantContext> implements NameAwareCodec {
    public final static int MAX_ALLELE_SIZE_BEFORE_WARNING = (int)Math.pow(2, 20);

    protected final static Logger log = Logger.getLogger(VCFCodec.class);
    protected final static int NUM_STANDARD_FIELDS = 8;  // INFO is the 8th column

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
    protected final String[] locParts = new String[6];

    // for performance we cache the hashmap of filter encodings for quick lookup
    protected HashMap<String,LinkedHashSet<String>> filterHash = new HashMap<String,LinkedHashSet<String>>();

    // we store a name to give to each of the variant contexts we emit
    protected String name = "Unknown";

    protected int lineNo = 0;

    protected Map<String, String> stringCache = new HashMap<String, String>();

    protected AbstractVCFCodec() {
        super(VariantContext.class);
    }

    /**
     * Creates a LazyParser for a LazyGenotypesContext to use to decode
     * our genotypes only when necessary.  We do this instead of eagarly
     * decoding the genotypes just to turn around and reencode in the frequent
     * case where we don't actually want to manipulate the genotypes
     */
    class LazyVCFGenotypesParser implements LazyGenotypesContext.LazyParser {
        final List<Allele> alleles;
        final String contig;
        final int start;

        LazyVCFGenotypesParser(final List<Allele> alleles, final String contig, final int start) {
            this.alleles = alleles;
            this.contig = contig;
            this.start = start;
        }

        @Override
        public LazyGenotypesContext.LazyData parse(final Object data) {
            //System.out.printf("Loading genotypes... %s:%d%n", contig, start);
            return createGenotypeMap((String) data, alleles, contig, start);
        }
    }

    /**
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    public abstract Object readHeader(LineReader reader);

    /**
     * create a genotype map
     *
     * @param str the string
     * @param alleles the list of alleles
     * @param chr chrom
     * @param pos position
     * @return a mapping of sample name to genotype object
     */
    public abstract LazyGenotypesContext.LazyData createGenotypeMap(String str, List<Allele> alleles, String chr, int pos);


    /**
     * parse the filter string, first checking to see if we already have parsed it in a previous attempt
     * @param filterString the string to parse
     * @return a set of the filters applied
     */
    protected abstract Set<String> parseFilters(String filterString);

    /**
     * create a VCF header from a set of header record lines
     *
     * @param headerStrings a list of strings that represent all the ## and # entries
     * @return a VCFHeader object
     */
    protected VCFHeader parseHeaderFromLines( final List<String> headerStrings, final VCFHeaderVersion version ) {
        Set<VCFHeaderLine> metaData = new TreeSet<VCFHeaderLine>();
        Set<String> sampleNames = new LinkedHashSet<String>();
        int contigCounter = 0;
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
                    sampleNames.add(strings[arrayIndex++]);

                if ( sawFormatTag && sampleNames.size() == 0 )
                    throw new UserException.MalformedVCFHeader("The FORMAT field was provided but there is no genotype/sample data");

            } else {
                if ( str.startsWith(VCFConstants.INFO_HEADER_START) ) {
                    final VCFInfoHeaderLine info = new VCFInfoHeaderLine(str.substring(7), version);
                    metaData.add(info);
                } else if ( str.startsWith(VCFConstants.FILTER_HEADER_START) ) {
                    final VCFFilterHeaderLine filter = new VCFFilterHeaderLine(str.substring(9), version);
                    metaData.add(filter);
                } else if ( str.startsWith(VCFConstants.FORMAT_HEADER_START) ) {
                    final VCFFormatHeaderLine format = new VCFFormatHeaderLine(str.substring(9), version);
                    metaData.add(format);
                } else if ( str.startsWith(VCFConstants.CONTIG_HEADER_START) ) {
                    final VCFContigHeaderLine contig = new VCFContigHeaderLine(str.substring(9), version, VCFConstants.CONTIG_HEADER_START.substring(2), contigCounter++);
                    metaData.add(contig);
                } else if ( str.startsWith(VCFConstants.ALT_HEADER_START) ) {
                    final VCFSimpleHeaderLine alt = new VCFSimpleHeaderLine(str.substring(6), version, VCFConstants.ALT_HEADER_START.substring(2), Arrays.asList("ID", "Description"));
                    metaData.add(alt);
                } else {
                    int equals = str.indexOf("=");
                    if ( equals != -1 )
                        metaData.add(new VCFHeaderLine(str.substring(2, equals), str.substring(equals+1)));
                }
            }
        }

        return new VCFHeader(metaData, sampleNames);
    }

    /**
     * the fast decode function
     * @param line the line of text for the record
     * @return a feature, (not guaranteed complete) that has the correct start and stop
     */
    public Feature decodeLoc(String line) {

        // the same line reader is not used for parsing the header and parsing lines, if we see a #, we've seen a header line
        if (line.startsWith(VCFHeader.HEADER_INDICATOR)) return null;

        // our header cannot be null, we need the genotype sample names and counts
        if (header == null) throw new ReviewedStingException("VCF Header cannot be null when decoding a record");

        final int nParts = ParsingUtils.split(line, locParts, VCFConstants.FIELD_SEPARATOR_CHAR, true);

        if ( nParts != 6 )
            throw new UserException.MalformedVCF("there aren't enough columns for line " + line, lineNo);

        // get our alleles (because the end position depends on them)
        final String ref = getCachedString(locParts[3].toUpperCase());
        final String alts = getCachedString(locParts[4].toUpperCase());
        final List<Allele> alleles = parseAlleles(ref, alts, lineNo);

        // find out our location
        int start = 0;
        try {
            start = Integer.valueOf(locParts[1]);
        } catch (Exception e) {
            generateException("the value in the POS field must be an integer but it was " + locParts[1], lineNo);
        }
        int stop = start;

        // ref alleles don't need to be single bases for monomorphic sites
        if ( alleles.size() == 1 ) {
            stop = start + alleles.get(0).length() - 1;
        }
        // we need to parse the INFO field to check for an END tag if it's a symbolic allele
        else if ( alleles.size() == 2 && alleles.get(1).isSymbolic() ) {
            final String[] extraParts = new String[4];
            final int nExtraParts = ParsingUtils.split(locParts[5], extraParts, VCFConstants.FIELD_SEPARATOR_CHAR, true);
            if ( nExtraParts < 3 )
                throw new UserException.MalformedVCF("there aren't enough columns for line " + line, lineNo);

            final Map<String, Object> attrs = parseInfo(extraParts[2]);
            try {
                stop = attrs.containsKey(VCFConstants.END_KEY) ? Integer.valueOf(attrs.get(VCFConstants.END_KEY).toString()) : start;
            } catch (Exception e) {
                throw new UserException.MalformedVCF("the END value in the INFO field is not valid for line " + line, lineNo);
            }
        }
        // handle multi-positional events
        else if ( !isSingleNucleotideEvent(alleles) ) {
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
    public VariantContext decode(String line) {
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
        VariantContextBuilder builder = new VariantContextBuilder();
        builder.source(getName());

        // increment the line count
        // TODO -- because of the way the engine utilizes Tribble, we can parse a line multiple times (especially when
        // TODO --   the first record is far along the contig) and the line counter can get out of sync
        lineNo++;

        // parse out the required fields
        final String chr = getCachedString(parts[0]);
        builder.chr(chr);
        int pos = Integer.valueOf(parts[1]);
        builder.start(pos);

        if ( parts[2].length() == 0 )
            generateException("The VCF specification requires a valid ID field");
        else if ( parts[2].equals(VCFConstants.EMPTY_ID_FIELD) )
            builder.noID();
        else
            builder.id(parts[2]);

        String ref = getCachedString(parts[3].toUpperCase());
        String alts = getCachedString(parts[4].toUpperCase());
        builder.log10PError(parseQual(parts[5]));
        builder.filters(parseFilters(getCachedString(parts[6])));
        final Map<String, Object> attrs = parseInfo(parts[7]);
        builder.attributes(attrs);

        // get our alleles, filters, and setup an attribute map
        List<Allele> alleles = parseAlleles(ref, alts, lineNo);

        // find out our current location, and clip the alleles down to their minimum length
        int stop = pos;
        // ref alleles don't need to be single bases for monomorphic sites
        if ( alleles.size() == 1 ) {
            stop = pos + alleles.get(0).length() - 1;
        }
        // we need to parse the INFO field to check for an END tag if it's a symbolic allele
        else if ( alleles.size() == 2 && alleles.get(1).isSymbolic() && attrs.containsKey(VCFConstants.END_KEY) ) {
            try {
                stop = Integer.valueOf(attrs.get(VCFConstants.END_KEY).toString());
            } catch (Exception e) {
                generateException("the END value in the INFO field is not valid");
            }
        } else if ( !isSingleNucleotideEvent(alleles) ) {
            ArrayList<Allele> newAlleles = new ArrayList<Allele>();
            stop = clipAlleles(pos, ref, alleles, newAlleles, lineNo);
            alleles = newAlleles;
        }
        builder.stop(stop);
        builder.alleles(alleles);

        // do we have genotyping data
        if (parts.length > NUM_STANDARD_FIELDS) {
            final LazyGenotypesContext.LazyParser lazyParser = new LazyVCFGenotypesParser(alleles, chr, pos);
            final int nGenotypes = header.getGenotypeSamples().size();
            LazyGenotypesContext lazy = new LazyGenotypesContext(lazyParser, parts[8], nGenotypes);

            // did we resort the sample names?  If so, we need to load the genotype data
            if ( !header.samplesWereAlreadySorted() )
                lazy.decode();

            builder.genotypesNoValidation(lazy);
        }

        VariantContext vc = null;
        try {
            builder.referenceBaseForIndel(ref.getBytes()[0]);
            vc = builder.make();
        } catch (Exception e) {
            generateException(e.getMessage());
        }

        return vc;
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
     * @return a mapping of keys to objects
     */
    private Map<String, Object> parseInfo(String infoField) {
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
        final int i;
        try {
            i = Integer.valueOf(index);
        } catch ( NumberFormatException e ) {
            throw new TribbleException.InternalCodecException("The following invalid GT allele index was encountered in the file: " + index);
        }
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
            return VariantContext.NO_LOG10_PERROR;

        Double val = Double.valueOf(qualString);

        // check to see if they encoded the missing qual score in VCF 3 style, with either the -1 or -1.0.  check for val < 0 to save some CPU cycles
        if ((val < 0) && (Math.abs(val - VCFConstants.MISSING_QUALITY_v3_DOUBLE) < VCFConstants.VCF_ENCODING_EPSILON))
            return VariantContext.NO_LOG10_PERROR;

        // scale and return the value
        return val / -10.0;
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

        if ( MAX_ALLELE_SIZE_BEFORE_WARNING != -1 && allele.length() > MAX_ALLELE_SIZE_BEFORE_WARNING )
            log.warn(String.format("Allele detected with length %d exceeding max size %d at approximately line %d, likely resulting in degraded VCF processing performance", allele.length(), MAX_ALLELE_SIZE_BEFORE_WARNING, lineNo));

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
     * return true if this is a symbolic allele (e.g. <SOMETAG>) or
     * structural variation breakend (with [ or ]), otherwise false
     * @param allele the allele to check
     * @return true if the allele is a symbolic allele, otherwise false
     */
    private static boolean isSymbolicAllele(String allele) {
        return (allele != null && allele.length() > 2 &&
                ((allele.startsWith("<") && allele.endsWith(">")) ||
                        (allele.contains("[") || allele.contains("]"))));
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

    public static boolean isSingleNucleotideEvent(List<Allele> alleles) {
        for ( Allele a : alleles ) {
            if ( a.length() != 1 )
                return false;
        }
        return true;
    }

    public static int computeForwardClipping(final List<Allele> unclippedAlleles, final byte ref0) {
        boolean clipping = true;
        int symbolicAlleleCount = 0;

        for ( Allele a : unclippedAlleles ) {
            if ( a.isSymbolic() ) {
                symbolicAlleleCount++;
                continue;
            }

            if ( a.length() < 1 || (a.getBases()[0] != ref0) ) {
                clipping = false;
                break;
            }
        }

        // don't clip if all alleles are symbolic
        return (clipping && symbolicAlleleCount != unclippedAlleles.size()) ? 1 : 0;
    }

    public static int computeReverseClipping(final List<Allele> unclippedAlleles, final byte[] ref, final int forwardClipping, final boolean allowFullClip, final int lineNo) {
        int clipping = 0;
        boolean stillClipping = true;

        while ( stillClipping ) {
            for ( final Allele a : unclippedAlleles ) {
                if ( a.isSymbolic() )
                    continue;

                // we need to ensure that we don't reverse clip out all of the bases from an allele because we then will have the wrong
                // position set for the VariantContext (although it's okay to forward clip it all out, because the position will be fine).
                if ( a.length() - clipping == 0 )
                    return clipping - (allowFullClip ? 0 : 1);

                if ( a.length() - clipping <= forwardClipping || a.length() - forwardClipping == 0 ) {
                    stillClipping = false;
                }
                else if ( ref.length == clipping ) {
                    if ( allowFullClip )
                        stillClipping = false;
                    else
                        generateException("bad alleles encountered", lineNo);
                }
                else if ( a.getBases()[a.length()-clipping-1] != ref[ref.length-clipping-1] ) {
                    stillClipping = false;
                }
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
    public static int clipAlleles(int position, String ref, List<Allele> unclippedAlleles, List<Allele> clippedAlleles, int lineNo) {

        int forwardClipping = computeForwardClipping(unclippedAlleles, (byte)ref.charAt(0));
        int reverseClipping = computeReverseClipping(unclippedAlleles, ref.getBytes(), forwardClipping, false, lineNo);

        if ( clippedAlleles != null ) {
            for ( Allele a : unclippedAlleles ) {
                if ( a.isSymbolic() ) {
                    clippedAlleles.add(a);
                } else {
                    final byte[] allele = Arrays.copyOfRange(a.getBases(), forwardClipping, a.getBases().length-reverseClipping);
                    if ( !Allele.acceptableAlleleBases(allele) )
                        generateException("Unparsable vcf record with bad allele [" + allele + "]", lineNo);
                    clippedAlleles.add(Allele.create(allele, a.isReference()));
                }
            }
        }

        // the new reference length
        int refLength = ref.length() - reverseClipping;

        return position+Math.max(refLength - 1,0);
    }

    public final static boolean canDecodeFile(final String potentialInput, final String MAGIC_HEADER_LINE) {
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
