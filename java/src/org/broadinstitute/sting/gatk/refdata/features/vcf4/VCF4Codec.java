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
public class VCF4Codec implements FeatureCodec {

    // we have to store the list of strings that make up the header until they're needed
    private VCFHeader header = null;

    // used to convert the index of the alternate allele in genotypes to a integer index
    private static int ZERO_CHAR = (byte)'0';

    // a mapping of the allele
    private static Map<String, List<Allele>> alleleMap = new HashMap<String, List<Allele>>(3);

    // cache the genotyope values
    private static String[] CachedGTValues = new String[100];

    // for performance testing purposes
    public static boolean validate = true;

    // a key optimization -- we need a per thread string parts array, so we don't allocate a big array over and over
    // todo: make this thread safe?
    private String[] parts = null;

    // for performance we cache the hashmap of filter encodings for quick lookup
    private HashMap<String,LinkedHashSet<String>> filterHash = new HashMap<String,LinkedHashSet<String>>();

    // a set of the genotype keys?
    private String[] genotypeKeyArray = new String[100];

    // a list of the info fields, filter fields, and format fields, for quick lookup to validate against
    ArrayList<String> infoFields = new ArrayList<String>();
    ArrayList<String> formatFields = new ArrayList<String>();
    ArrayList<String> filterFields = new ArrayList<String>();

    // do we want to validate the info, format, and filter fields
    private final boolean validateFromHeader = true;

    // we store a name to give to each of the variant contexts we emit
    private String name = "Unkown";

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
    private int createHeader(List<String> headerStrings, String line) {
        headerStrings.add(line);
        header = VCFReaderUtils.createHeader(headerStrings, VCFHeaderVersion.VCF4_0);

        // load the parsing fields
        Set<VCFHeaderLine> headerLines = header.getMetaData();

        // setup our look-up lists for validation
        for (VCFHeaderLine hl : headerLines) {
            if (hl.getClass() == VCFFilterHeaderLine.class)
                this.filterFields.add(((VCFFilterHeaderLine)hl).getmName());
            if (hl.getClass() == VCFFormatHeaderLine.class)
                                        this.formatFields.add(((VCFFormatHeaderLine)hl).getmName());
            if (hl.getClass() == VCFInfoHeaderLine.class)
                                        this.infoFields.add(((VCFInfoHeaderLine)hl).getmName());
        }
        // sort the lists so we can binary search them later on
        Collections.sort(filterFields);
        Collections.sort(formatFields);
        Collections.sort(infoFields);

        return headerStrings.size();
    }

    /**
     * the fast decode function
     * @param line the line of text for the record
     * @return a feature, (not guaranteed complete) that has the correct start and stop
     */
    public Feature decodeLoc(String line) {
        return decode(line);
    }

    /**
     * decode the line into a feature (VariantContext)
     * @param line the line
     * @return a VariantContext
     */
    public Feature decode(String line) {
        if (parts == null)
            parts = new String[header.getColumnCount()];

        int nParts = ParsingUtils.split(line, parts, '\t');

        // our header cannot be null, we need the genotype sample names and counts
        if (header == null) throw new IllegalStateException("VCF Header cannot be null");

        // check to make sure the split resulted in the correct number of fields (8 + (1 + genotytpe counts if it has genotypes)
        if (nParts != header.getColumnCount())
            throw new IllegalArgumentException("we expected " + header.getColumnCount() + " columns and we got " + nParts + " for line " + line);


        return parseVCFLine(parts);
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
    private static List<Allele> parseGenotypeAlleles(String GT, List<Allele> alleles, Map<String, List<Allele>> cache) {
        // this should cache results [since they are immutable] and return a single object for each genotype
        if ( GT.length() != 3 ) throw new StingException("Unreasonable number of alleles"); // 0/1 => barf on 10/0
        List<Allele> GTAlleles = cache.get(GT);
        if ( GTAlleles == null ) {
            GTAlleles = Arrays.asList(oneAllele(GT.charAt(0), alleles), oneAllele(GT.charAt(2), alleles));
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
                    value = field.substring(eqI+1, field.length()); // todo -- needs to convert to int, double, etc
                    //System.out.printf("%s %s%n", key, value);
                } else {
                    key = field;
                    value = 1;
                }

                attributes.put(key, value);
            }
        }
        // validate the fields
        validateFields(attributes.keySet(),infoFields);

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
                    throw new StingException("Unable to find field descibing attribute " + attr);
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
        checkAllele(ref, true);
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
     */
    private static void checkAllele(String allele,boolean isRef) {
        if ( ! Allele.acceptableAlleleBases(allele,isRef) ) {
            throw new StingException("Unparsable vcf record with allele " + allele);
        }
    }

    /**
     * parse a single allele, given the allele list
     * @param alleles the alleles available
     * @param alt the allele to parse
     * @param isRef are we the reference allele?
     */
    private void parseSingleAllele(List<Allele> alleles, String alt, boolean isRef) {
        checkAllele(alt,isRef);

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
    private VariantContext parseVCFLine(String[] parts) {
        // parse out the required fields
        String contig = parts[0];
        long pos = Long.valueOf(parts[1]);
        String id = parts[2];
        String ref = parts[3].toUpperCase();
        String alts = parts[4].toUpperCase();
        Double qual = parseQual(parts[5]);
        String filter = parts[6];
        String info = parts[7];


        List<Allele> alleles = parseAlleles(ref, alts);
        Set<String> filters = parseFilters(filter);
        Map<String, Object> attributes = parseInfo(info, id);

        // find out our current location, and clip the alleles down to their minimum length
        Pair<GenomeLoc, List<Allele>> locAndAlleles = (ref.length() > 1) ?
                clipAlleles(contig, pos, ref, alleles) :
                new Pair<GenomeLoc, List<Allele>>(GenomeLocParser.createGenomeLoc(contig, pos), alleles);

        // a map to store our genotypes
        Map<String, Genotype> genotypes = null;

        // do we have genotyping data
        if (parts.length > 8) {
            String GT = parts[8];
            int genotypesStart = 9;
            // parse genotypes
            int nGTKeys = ParsingUtils.split(GT, genotypeKeyArray, ':');
            genotypes = new HashMap<String, Genotype>(Math.max(parts.length - genotypesStart, 1));
            Iterator<String> iter = header.getGenotypeSamples().iterator();

            alleleMap.clear();
            for (int genotypeOffset = genotypesStart; genotypeOffset < parts.length; genotypeOffset++) {
                String sample = parts[genotypeOffset];
                String[] GTValues = CachedGTValues;
                ParsingUtils.split(sample, GTValues, ':');
                List<Allele> genotypeAlleles = parseGenotypeAlleles(GTValues[0], locAndAlleles.second, alleleMap);
                double GTQual = VariantContext.NO_NEG_LOG_10PERROR;
                Set<String> genotypeFilters = null;

                // todo -- the parsing of attributes could be made lazy for performance
                Map<String, String> gtAttributes = null;
                if (nGTKeys > 1) {
                    gtAttributes = new HashMap<String, String>(nGTKeys - 1);
                    for (int i = 1; i < nGTKeys; i++) {
                        if (genotypeKeyArray[i].equals("GQ")) {
                            GTQual = parseQual(GTValues[i]);
                        }
                        if (genotypeKeyArray[i].equals("FL")) { // deal with genotype filters here
                            genotypeFilters.addAll(parseFilters(GTValues[i]));
                        } else {
                            gtAttributes.put(genotypeKeyArray[i], GTValues[i]);
                        }
                    }
                    // validate the format fields
                    validateFields(gtAttributes.keySet(), formatFields);
                }

                boolean phased = genotypeKeyArray[0].charAt(1) == '|';

                Genotype g = new Genotype(iter.next(), genotypeAlleles, GTQual, genotypeFilters, gtAttributes, phased);
                genotypes.put(g.getSampleName(), g);
            }
        }
        // todo -- we need access to our track name to name the variant context
        return new VariantContext(name, locAndAlleles.first, locAndAlleles.second, genotypes, qual, filters, attributes);
    }

    /**
     * clip the alleles, based on the reference
     * @param unclippedAlleles the list of alleles
     * @return a list of alleles, clipped to the reference
     */
    static Pair<GenomeLoc,List<Allele>> clipAlleles(String contig, long position, String ref, List<Allele> unclippedAlleles) {
        List<Allele> newAlleleList = new ArrayList<Allele>();

        // find the preceeding string common to all alleles and the reference
        int forwardClipped = 0;
        boolean clipping = true;
        while (clipping) {
            for (Allele a : unclippedAlleles)
                if (forwardClipped > ref.length() - 1)
                    clipping = false;
                else if (a.length() <= forwardClipped || (a.getBases()[forwardClipped] != ref.getBytes()[forwardClipped])) {
                    clipping = false;
                }
            if (clipping) forwardClipped++;
        }

        int reverseClipped = 0;
        clipping = true;
        while (clipping) {
            for (Allele a : unclippedAlleles)
                if (a.length() - reverseClipped < 0 || a.length() - forwardClipped == 0)
                    clipping = false;
                else if (a.getBases()[a.length()-reverseClipped-1] != ref.getBytes()[ref.length()-reverseClipped-1])
                    clipping = false;
            if (clipping) reverseClipped++;
        }

        // check to see if we're about to clip all the bases from the reference, if so back off the front clip a base
        if (forwardClipped + reverseClipped >= ref.length())
            forwardClipped--;

        for (Allele a : unclippedAlleles)
            newAlleleList.add(Allele.create(Arrays.copyOfRange(a.getBases(),forwardClipped,a.getBases().length-reverseClipped),a.isReference()));
        return new Pair<GenomeLoc,List<Allele>>(GenomeLocParser.createGenomeLoc(contig,position+forwardClipped,(position+ref.length()-reverseClipped-1)),
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
